"""Modulation Transfer Function Module for slanted edge figures

The optical transfer function (OTF) of an optical system specifies 
how different spatial frequencies are handled by the system. 
It is used by optical engineers to describe how the optics project 
light from the object or scene onto detector or simply the next item 
in the optical transmission chain. 
A variant, the modulation transfer function (MTF), neglects phase effects, 
but is equivalent to the OTF in many situations.

look also: https://en.wikipedia.org/wiki/Optical_transfer_function

In normal operation MTF is fourier transform of Point Spread Function.
But for slanted edge targeted MTF calculations we use LSF instead of PSF.
This is because we assume LSF as cross-section of PSF along slant.
LSF can be optained using detivative of Edge Spread Function (ESF).

For line scan operations it is better to calculate along track and 
across track MTF's.

Also for non symetrical pixel pitch/sizes it will also be convinient
to calculate horizontal and vertical slants.

Note that, the prefered slant tilt is between 2 degree and 10 degree.

Dependencies
    Pillow numpy scipy matpilotlib opencv-python
"""
import matplotlib.pyplot as plt
import pylab as pylab
import numpy as np
import cv2 as cv2
import math as math

from PIL import Image, ImageOps
from scipy import interpolate
from scipy.fft import fft
from enum import Enum
from dataclasses import dataclass

@dataclass
class cSet:
    """
    Base class to store index value pairs
    ...
    Attributes
    x : np.ndarray
        numpy array to store indexes
    y : np.ndarray
        numpy array to store values
    """    
    x: np.ndarray
    y: np.ndarray

@dataclass
class cESF:
    """
    Base class to store Edge Spread Function (ESF) analysis output
    ...
    Attributes
    rawESF : cSet
        ESF optained from image
    interpESF : cSet
        Interpolated ESF of rawESF
    threshold :float
        threshold used to get ESF
    width : float
        pixel transition width
    angle : float
        slant angle in degrees
    edgePoly : np.ndarray
        polynomial of edge slant
    """    
    rawESF: cSet
    interpESF: cSet
    threshold:float
    width: float
    angle: float
    edgePoly: np.ndarray

@dataclass
class cMTF:
    """
    Base class to Modulation Transfer Function (MTF) values
    ...
    Attributes
    x : np.ndarray
        numpy array to store indexes
    y : np.ndarray
        numpy array to store values
    mtfAtNyquist : float
        MTF value at Nyquist Frequency
    width : float
        Pixel transition
    """    
    x: np.ndarray
    y: np.ndarray
    mtfAtNyquist: float
    width: float

class Verbosity(Enum):
    """
    Output types of module methods
    ...
    Enumerations
    NONE : no output
    BRIEF : brif text output
    DETAIL : graphical output
    """    
    NONE = 0
    BRIEF = 1
    DETAIL = 2

class Helper:
    @staticmethod
    def LoadImage(filename):
        """
        Load image from given path
        ...
        Parameters
        filename : str
            a fully quialified file name of the image
        ...
        Returns
        PIL.Image
            Grayscale image data
        """    
        img = Image.open(filename)
        if img.mode in {'I;16','I;16L','I;16B','I;16N'}:
            gsimg = img
        else:
            gsimg = img.convert('L')
        return gsimg

    @staticmethod
    def LoadImageAsArray(filename):
        """
        Load and convert image from given path to numpy array
        ...
        Parameters
        filename : str
            a fully quialified file name of the image
        ...
        Returns
        np.darray
            Grayscale image data as numpy array with values between 0.0 and 1.0
        """    
        img = Helper.LoadImage(filename)
        if img.mode in {'I;16','I;16L','I;16B','I;16N'}:
            arr = np.asarray(img, dtype=np.double)/65535
        else:
            arr = np.asarray(img, dtype=np.double)/255
        return arr

    @staticmethod
    def ImageToArray(img):
        """
        Convert PIL image to numpy array
        ...
        Parameters
        img : PIL.Image
            Source image
        ...
        Returns
        np.darray
            Grayscale image data as numpy array with values between 0.0 and 1.0
        """    
        if img.mode in {'I;16','I;16L','I;16B','I;16N'}:
            arr = np.asarray(img, dtype=np.double)/65535
        else:
            arr = np.asarray(img, dtype=np.double)/255
        return arr

    @staticmethod
    def ArrayToImage(imgArr):
        """
        Convert numpy array to PIL image
        ...
        Parameters
        imgArr : np.darray
            Image data as numpy array having values between 0.0 and 1.0
        ...
        Returns
        PIL.Image
            Grayscale PIL image
        """    
        img = Image.fromarray(imgArr*255, mode='L')
        return img

    @staticmethod
    def CorrectImageOrientation(imgArr):
        """
        Rotate, transpose, flip image to get correct orientation for analysis
        Module assumes a "dark side up" horizontal slanted image
        ...
        Parameters
        imgArr : np.darray
            Image data as numpy array having values between 0.0 and 1.0
        ...
        Returns
        np.darray
            Orientation corrected image array
        """    
        tl = np.average(imgArr[0:2,0:2])
        tr = np.average(imgArr[0:2,-3:-1])
        bl = np.average(imgArr[-3:-1,0:2])
        br = np.average(imgArr[-3:-1,-3:-1])
        edges = [tl, tr, bl, br]
        edgeIndexes = np.argsort(edges)
        if (edgeIndexes[0] + edgeIndexes[1]) == 1:
            pass
        elif (edgeIndexes[0] + edgeIndexes[1]) == 5:
            imgArr = np.flip(imgArr, axis=0)
        elif (edgeIndexes[0] + edgeIndexes[1]) == 2:
            imgArr = np.transpose(imgArr)
        elif (edgeIndexes[0] + edgeIndexes[1]) == 4:
            imgArr = np.flip(np.transpose(imgArr), axis=0)
        
        return imgArr
  
class MTF:
    @staticmethod
    def SafeCrop(values, distances, head, tail):
        """
        Safely crop a index-value array from head to tail
        Note that, method does not crop operation uses values of distance
          not index itself
        ...
        Parameters
        value : np.darray
            Value array
        distances : np.darray
            Index array
        head : float
            Desired crop start index value
        tail : float
            Desired crop end index value
        ...
        Returns
        cSet
            crop result as cSet
        """
        isIncrementing = True
        if distances[0] > distances[-1]:
            isIncrementing = False
            distances = -distances
            dummy = -tail
            tail = -head
            head = dummy


        hindex = (np.where(distances < head)[0])
        tindex = (np.where(distances > tail)[0])

        if hindex.size < 2:
            h = 0
        else:
            h = np.amax(hindex)

        if tindex.size == 0:
            t = distances.size
        else:
            t = np.amin(tindex)

        if isIncrementing == False:
            distances = -distances

        return cSet(distances[h:t], values[h:t])

    @staticmethod
    def GetEdgeSpreadFunction(imgArr, edgePoly, verbose=Verbosity.NONE):
        """
        Calculate Edge Spred Function (ESF)
        
        ESF is simply distance map of every pixel to the given edge polynomial
        https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line
        ...
        Parameters
        imgArr : np.darray
            Image data as numpy array having values between 0.0 and 1.0
        edgePoly : np.darray
            Edge polynomial parameters
        verbose : Verbosity
            Output verbosity level
        ...
        Returns
        cSet
            ESF as value index set
        """
        Y = imgArr.shape[0]
        X = imgArr.shape[1]

        values = np.reshape(imgArr, X*Y)

        distance = np.zeros((Y,X))
        column = np.arange(0,X)+0.5
        for y in range(Y):
            distance[y,:] = (edgePoly[0]*column - (y+0.5) + edgePoly[1]) / np.sqrt(edgePoly[0]*edgePoly[0] + 1)

        distances = np.reshape(distance, X*Y)
        indexes = np.argsort(distances)
        
        sign = 1
        if np.average(values[indexes[:10]]) > np.average(values[indexes[-10:]]):
            sign = -1

        values = values[indexes]
        distances = sign*distances[indexes]     
        
        if (distances[0] > distances[-1]):
            distances = np.flip(distances)
            values = np.flip(values)

        if (verbose == Verbosity.BRIEF):
            print("Raw ESF [done] (Distance from {0:2.2f} to {1:2.2f})".format(sign*distances[0], sign*distances[-1]))

        elif (verbose == Verbosity.DETAIL):
            x = [0, np.size(imgArr,1)-1]
            y = np.polyval(edgePoly, x)

            fig = pylab.gcf()
            fig.canvas.manager.set_window_title('Raw ESF')
            (ax1, ax2) = plt.subplots(2)
            ax1.imshow(imgArr, cmap='gray', vmin=0.0, vmax=1.0)
            ax1.plot(x, y, color='red')
            ax2.plot(distances, values)
            plt.show()
            plt.show(block=False)

        return cSet(distances, values)

    @staticmethod
    def GetEdgeSpreadFunctionCrop(imgArr, verbose=Verbosity.NONE):
        """
        Calculate and crop Edge Spread Function (ESF). 
        Crop occures around center of transition
        ...
        Parameters
        imgArr : np.darray
            Image data as numpy array having values between 0.0 and 1.0
        verbose : Verbosity
            Output verbosity level
        ...
        Returns
        cESF
            Extended ESF information
        """
        imgArr = Helper.CorrectImageOrientation(imgArr)
        edgeImg = cv2.Canny(np.uint8(imgArr*255), 40, 90, L2gradient=True)

        line = np.argwhere(edgeImg == 255)
        edgePoly = np.polyfit(line[:,1],line[:,0],1)
        angle = math.degrees(math.atan(-edgePoly[0]))

        finalEdgePoly = edgePoly.copy()
        if angle > 0:
            imgArr = np.flip(imgArr, axis=1)
            finalEdgePoly[1] = np.polyval(edgePoly,np.size(imgArr,1)-1)
            finalEdgePoly[0] = -edgePoly[0]

        esf = MTF.GetEdgeSpreadFunction(imgArr, finalEdgePoly, Verbosity.NONE)

        esfValues = esf.y
        esfDistances = esf.x

        maximum = np.amax(esfValues)
        minimum = np.amin(esfValues)

        threshold = (maximum - minimum) * 0.1

        head = np.amax(esfDistances[(np.where(esfValues < minimum + threshold))[0]])
        tail = np.amin(esfDistances[(np.where(esfValues > maximum - threshold))[0]])

        width = abs(head-tail)

        esfRaw = MTF.SafeCrop(esfValues, esfDistances, head - 1.2*width, tail + 1.2*width)

        qs = np.linspace(0,1,20)[1:-1]
        knots = np.quantile(esfRaw.x, qs)
        tck = interpolate.splrep(esfRaw.x, esfRaw.y, t=knots, k=3)
        ysmooth = interpolate.splev(esfRaw.x, tck)
        
        InterpDistances = np.linspace(esfRaw.x[0], esfRaw.x[-1], 500)
        InterpValues = np.interp(InterpDistances, esfRaw.x, ysmooth)
        
        esfInterp = cSet(InterpDistances, InterpValues)

        if (verbose == Verbosity.BRIEF):
            print("ESF Crop [done] (Distance from {0:2.2f} to {1:2.2f})".format(esfRaw.x[0], esfRaw.x[-1]))

        elif (verbose == Verbosity.DETAIL):
            x = [0, np.size(imgArr,1)-1]
            y = np.polyval(finalEdgePoly, x)

            fig = pylab.gcf()
            fig.canvas.manager.set_window_title('ESF Crop')
            (ax1, ax2) = plt.subplots(2)
            ax1.imshow(imgArr, cmap='gray', vmin=0.0, vmax=1.0)
            ax1.plot(x, y, color='red')
            ax2.plot(esfRaw.x, esfRaw.y,InterpDistances,InterpValues)
            plt.show(block=False)
            plt.show()

        return cESF(esfRaw, esfInterp, threshold, width, angle, edgePoly)

    @staticmethod
    def SimplifyEdgeSpreadFunction(esf, verbose=Verbosity.NONE):
        """
        Remove dublicate distance occurances of Edge Spread Function. 
        ...
        Parameters
        esf : cSet
            ESF index-value set
        verbose : Verbosity
            Output verbosity level
        ...
        Returns
        cSet
            Shrinked ESF index-value set
        """
        res = np.unique(esf.x, return_index=True, return_counts=True)

        indexes = res[1]
        counts = res[2]
        sz = np.size(res[0])
        
        distances = esf.x[indexes]
        values = np.zeros(sz, dtype=np.float)
        
        for x in range(sz):
            values[x] = np.sum(esf.y[indexes[x]:indexes[x]+counts[x]])/counts[x]

        if (verbose == Verbosity.BRIEF):
            print("ESF Simplification [done] (Size from {0:d} to {1:d})".format(np.size(esf.x), np.size(distances)))

        elif (verbose == Verbosity.DETAIL):
            fig = pylab.gcf()
            fig.canvas.manager.set_window_title("ESF Simplification (Size from {0:d} to {1:d})".format(np.size(esf.x), np.size(distances)))
            (ax1, ax2) = plt.subplots(2)
            ax1.plot(esf.x, esf.y)
            ax2.plot(distances, values)
            plt.show(block=False)
            plt.show()

        return cSet(distances, values)

    @staticmethod
    def GetLineSpreadFunction(esf, normalize=True, verbose=Verbosity.NONE):
        """
        Calculate Line Spread Function (LSF) from ESF
        
        LSF is simply the derivative of ESF.
        For a better result, an interpolated ESF can be used instead of raw ESF
        ...
        Parameters
        esf : cSet
            ESF index-value set
        normalize : bool
            Calculated normalized LSF
        verbose : Verbosity
            Output verbosity level
        ...
        Returns
        cSet
            LSF index-value set
        """
        lsfDividend = np.diff(esf.y)
        lsfDivisor  = np.diff(esf.x)

        lsfValues = np.divide(lsfDividend, lsfDivisor)
        lsfDistances = esf.x[0:-1]

        if normalize:
            lsfValues = lsfValues / (max(lsfValues))

        if (verbose == Verbosity.BRIEF):
            print("MTF [done]")

        elif (verbose == Verbosity.DETAIL):
            fig = pylab.gcf()
            fig.canvas.manager.set_window_title("LSF")
            (ax1) = plt.subplots(1)
            ax1.plot(lsfDistances, lsfValues)
            plt.show(block=False)
            plt.show()

        return cSet(lsfDistances, lsfValues)

    @staticmethod
    def GetMTF(lsf, verbose=Verbosity.NONE):
        """
        Calculate Modulation Transform Function (MTF) from LSF
        ...
        Parameters
        lsf : cSet
            LSF index-value set
        verbose : Verbosity
            Output verbosity level
        ...
        Returns
        cMTF
            MTF value set except pixel transition width
        """
        N = np.size(lsf.x)
        px = N/(lsf.x[-1]-lsf.x[0])
        values = 1/np.sum(lsf.y)*abs(fft(lsf.y))
        distances = np.arange(0,N)/N*px

        interpDistances = np.linspace(0,1,50)
        interp = interpolate.interp1d(distances, values, kind='cubic')
        interpValues = interp(interpDistances)
        valueAtNyquist = interpValues[25]*100

        if (verbose == Verbosity.BRIEF):
            print("MTF [done]")

        elif (verbose == Verbosity.DETAIL):
            fig = pylab.gcf()
            fig.canvas.manager.set_window_title("MTF ({0:2.2f}% at Nyquist)".format(valueAtNyquist))
            (ax1) = plt.subplots(1)
            ax1.plot(interpDistances, interpValues)
            #ax1.plot( values)
            plt.show(block=False)
            plt.show()
        
        return cMTF(interpDistances, interpValues, valueAtNyquist, -1.0)

    @staticmethod
    def CalculateMtf(imgArr, verbose=Verbosity.NONE):
        """
        Calculate Modulation Transform Function (MTF) of an image array
        ...
        Parameters
        imgArr : np.darray
            Image data as numpy array having values between 0.0 and 1.0
        verbose : Verbosity
            Output verbosity level
        ...
        Returns
        cMTF
            MTF value set 
        """
        imgArr = Helper.CorrectImageOrientation(imgArr)
        esf = MTF.GetEdgeSpreadFunctionCrop(imgArr, Verbosity.NONE)
        lsf = MTF.GetLineSpreadFunction(esf.interpESF, True, Verbosity.NONE)
        mtf = MTF.GetMTF(lsf, Verbosity.NONE)

        if (verbose == Verbosity.BRIEF):
            print("MTF at Nyquist:{0:0.2f}%, Transition Width:{1:0.2f}".format(mtf.mtfAtNyquist, esf.width))

        elif (verbose == Verbosity.DETAIL):
            x = [0, np.size(imgArr,1)-1]
            y = np.polyval(esf.edgePoly, x)

            fig = pylab.gcf()
            fig.canvas.manager.set_window_title('MTF Analysis')
            gs = fig.add_gridspec(3,2)
            ax1 = fig.add_subplot(gs[0, 0])
            ax2 = fig.add_subplot(gs[1, 0])
            ax3 = fig.add_subplot(gs[2, 0])
            ax4 = fig.add_subplot(gs[:, 1])

            ax1.imshow(imgArr, cmap='gray', vmin=0.0, vmax=1.0)
            ax1.plot(x, y, color='red')
            ax1.axis('off')
            ax2.plot(esf.rawESF.x, esf.rawESF.y,
                     esf.interpESF.x, esf.interpESF.y)
            top = np.max(esf.rawESF.y)-esf.threshold
            bot = np.min(esf.rawESF.y)+esf.threshold
            ax2.plot([esf.rawESF.x[0], esf.rawESF.x[-1]], [top, top], color='red')
            ax2.plot([esf.rawESF.x[0], esf.rawESF.x[-1]], [bot, bot], color='red')
            ax2.xaxis.set_visible(False)
            ax2.yaxis.set_visible(False)
            ax3.plot(lsf.x, lsf.y)
            ax3.xaxis.set_visible(False)
            ax3.yaxis.set_visible(False)
            ax4.plot(mtf.x, mtf.y)
            ax4.set_title("MTF at Nyquist:{0:0.2f}%\nTransition Width:{1:0.2f}".format(mtf.mtfAtNyquist, esf.width))
            ax4.grid(True)

            plt.show(block=False)
            plt.show()

        return cMTF(x, y, mtf.mtfAtNyquist, esf.width)
