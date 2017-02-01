import pyfits
from numpy import *
import sys
from DavidsNM import save_img

f = pyfits.open(sys.argv[1])
dat = f[0].data
f.close()

f = pyfits.open(sys.argv[2])
dat2 = f[0].data
f.close()

save_img(dat/dat2, "div.fits")
