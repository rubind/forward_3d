from numpy import *
import pyfits
import sys
from scipy import fftpack as ft
from DavidsNM import save_img

im = sys.argv[1]
wv = float(sys.argv[2])

f = pyfits.open(im)
resid = f[0].data
f.close()

f = pyfits.open("../../psfs/IFC_PSF_%.2f_micron_50_mas.fits" % wv)
psf = f[0].data
f.close()

resid = array(resid, dtype=float64)
psf = array(psf, dtype=float64)


conv = np.real(ft.ifft2(ft.fft2(resid) * ft.fft2(psf))

newim = im.replace(".fits", "_conv.fits")
assert newim != im

save_img(conv, newim)
