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

try:
    f = pyfits.open("../../psfs/IFC_PSF_%.2f_micron_50+150_mas.fits" % wv)
except:
    f = pyfits.open("../../../psfs/IFC_PSF_%.2f_micron_50+150_mas.fits" % wv)

psf = f[0].data
f.close()

resid = array(resid, dtype=float64)
psf = array(psf, dtype=float64)

print resid.shape

conv = resid*0

for i in range(conv.shape[0]):
    for j in range(conv.shape[1]):
        theshape = conv.shape[2:]
        

        psftmp = psf[psf.shape[0]/2-theshape[0]/2: psf.shape[0]/2+theshape[0]/2,
                     psf.shape[1]/2-theshape[1]/2: psf.shape[1]/2+theshape[1]/2]

        print psftmp.shape

        adjust = psftmp *0
        adjust[psftmp.shape[0]/2, psftmp.shape[1]/2] = 1.

        conv[i, j] = real(ft.ifft2(ft.fft2(resid[i, j]) * ft.fft2(psftmp) * ft.fft2(adjust)))

newim = im.replace(".fits", "_conv.fits")
assert newim != im

save_img(conv, newim)
