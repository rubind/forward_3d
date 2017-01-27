import pyfits
from numpy import *
import sys

newim = sys.argv[1].replace(".fits", "_ds9.fits")
assert newim != sys.argv[1]

f = pyfits.open(sys.argv[1])
dat = f[0].data


print dat.shape
dat = transpose(dat, [0,3,1,2])
print dat.shape

f[0].data = dat
f.writeto(newim, clobber = "True")
f.close()

