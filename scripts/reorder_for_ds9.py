import pyfits
from numpy import *
import sys

for i in range(1, len(sys.argv)):
    newim = sys.argv[i].replace(".fits", "_ds9.fits")
    assert newim != sys.argv[i]
    
    f = pyfits.open(sys.argv[i])
    dat = f[0].data
    
    
    print dat.shape
    dat = transpose(dat, [0,3,1,2])
    print dat.shape
    
    f[0].data = dat
    f.writeto(newim, clobber = "True")
    f.close()

