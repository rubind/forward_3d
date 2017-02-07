import pyfits
from numpy import *
import sys
from DavidsNM import save_img

if " ".join(sys.argv).count(".fits"):

    f = pyfits.open(sys.argv[1])
    dat = f[0].data
    f.close()

    f = pyfits.open(sys.argv[2])
    dat2 = f[0].data
    f.close()

    save_img(dat/dat2, "div.fits")

else:
    
    for dr in sys.argv[1:]:
        f = pyfits.open(dr + "/resid_ds9.fits")
        dat = f[0].data
        f.close()
        
        f = pyfits.open(dr + "/true_scene_ds9.fits")
        dat2 = f[0].data
        f.close()

        save_img(dat/dat2, dr + "/div.fits")
        


