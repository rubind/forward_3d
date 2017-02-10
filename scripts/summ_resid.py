import pyfits
from numpy import *
import sys

if len(sys.argv) == 1:
    wds = ["."]
else:
    wds = sys.argv[1:]

for wd in wds:
    f = pyfits.open(wd + "/resid.fits")
    dat = f[0].data
    f.close()
    
    f = pyfits.open(wd + "/true_scene.fits")
    dat2 = f[0].data
    f.close()
    
    frac_resid = dat/dat2
    inds = where(1  - isnan(frac_resid))
    frac_resid = frac_resid[inds]
    
    print std(frac_resid)
    print 1.4826*median(abs(frac_resid - median(frac_resid)))






    for dith in range(dat.shape[0]):
        
        frac_resid = abs(dat[dith])/(abs(dat2[dith]).max())
        
        print "max", dith, abs(dat2[dith]).max()
        print wd, "max", dith, abs(dat[dith]).max(), "=", abs(dat[dith]).max()/abs(dat2[dith]).max()
        
        inds = where(1  - isnan(frac_resid))
        frac_resid = frac_resid[inds]
        
        print dith, std(frac_resid)
        print dith, 1.4826*median(abs(frac_resid - median(frac_resid)))
