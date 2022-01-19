import mtf as mtf

imgArr = mtf.Helper.LoadImageAsArray('slant.png')
res = mtf.MTF.CalculateMtf(imgArr, verbose=mtf.Verbosity.DETAIL)


