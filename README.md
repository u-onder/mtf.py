<div id="top"></div>
<!-- PROJECT LOGO -->
<br />
<div align="center">
  <h2 align="center">MTF.py</h2>
  <h3 align="center">Modulation Transfer Function Calculation for Slanted Edge Images</h3>
</div>

<!-- ABOUT -->
## About

The optical transfer function (OTF) of an optical system specifies 
how different spatial frequencies are handled by the system. 
It is used by optical engineers to describe how the optics project 
light from the object or scene onto the detector or simply the next item 
in the optical transmission chain. 
A variant, the modulation transfer function (MTF), neglects phase effects, 
but is equivalent to the OTF in many situations.

In normal operation, MTF is the Fourier transform of the Point Spread Function (PSF).
But for slanted edge-targeted MTF calculations, we use Line Spread Function (LSF)
instead of PSF. This is because we assume LSF as a cross-section of PSF along the slant.
LSF can be obtained using the derivative of Edge Spread Function (ESF).

Since LSF is the cross-section of PSF along the slant edge, a better approach will be 
calculating MTF for both horizontally and vertically tilted images.

For line-scan operations, it is better to calculate along the track and 
across the track MTF's.

Also for non-symmetrical pixel pitch/sizes sensors, it will also be convenient
to calculate horizontal and vertical slants.

Note that the preferred slant tilt is between 2 degrees and 10 degrees.

<p align="right">(<a href="#top">back to top</a>)</p>

<!-- USAGE -->
## Usage
import mtf as mtf

imgArr = mtf.Helper.LoadImageAsArray('slant.png')
res = mtf.MTF.CalculateMtf(imgArr, verbose=mtf.Verbosity.DETAIL)

<p align="right">(<a href="#top">back to top</a>)</p>

<!-- DEPENDENCIES -->
## Dependencies

This module requires the following external modules
    Pillow, numpy, scipy, matplotlib, opencv-python

<p align="right">(<a href="#top">back to top</a>)</p>

<!-- LICENSE -->
## License

Distributed under the MIT License. See `LICENSE.txt` for more information.

<p align="right">(<a href="#top">back to top</a>)</p>

<!-- CONTACT -->
## Contact

Ufuk Onder - ufuk.onder@gmail.com

Project Link: (https://github.com/u-onder/mtf.py)

<p align="right">(<a href="#top">back to top</a>)</p>

<!-- ACKNOWLEDGMENTS -->
## Acknowledgments

* [Optical Transfer Function](https://en.wikipedia.org/wiki/Optical_transfer_function)
* [How to Measure Modulation Transfer Function](https://harvestimaging.com/blog/?p=1328)
* [Comparison of MTF Measurements Using Edge Method](https://hal.archives-ouvertes.fr/hal-02055611/document)

<p align="right">(<a href="#top">back to top</a>)</p>
