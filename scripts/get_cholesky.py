from numpy import *
import pyfits
from DavidsNM import save_img
import time

print time.asctime()
f = pyfits.open("jacobian.fits")
j = f[0].data
f.close()

print j.shape
print time.asctime()

cmat = dot(transpose(j), j)
j = 0

print cmat.shape
print time.asctime()

Lmat  = linalg.cholesky(cmat)

save_img(Lmat, "Lmat.fits")
print time.asctime()
print "Done!"
