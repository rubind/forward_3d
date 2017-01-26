import pyfits
from DavidsNM import save_img
from numpy import *

f = pyfits.open("jacobian.fits")
j = f[0].data
f.close()

f = pyfits.open("true_node_vals.fits")
true_node_vals = f[0].data
f.close()

f = pyfits.open("true_scene.fits")
true_scene = f[0].data
f.close()

true_scene_reshape = reshape(true_scene, len(j))

print j.shape, true_scene.shape

solved = linalg.lstsq(j, true_scene_reshape)[0]

est_model = dot(j, solved)

save_img(reshape(est_model, true_scene.shape), "resid.fits")

print abs(solved - true_node_vals).max()
