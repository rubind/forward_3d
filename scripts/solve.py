import pyfits
from DavidsNM import save_img
from numpy import *
import sys

regularization = float(sys.argv[1])

try:
    f = pyfits.open("jacobian.fits")
except:
    f = pyfits.open("../jacobian.fits")

j = f[0].data
f.close()

try:
    f = pyfits.open("true_node_vals.fits")
except:
    f = pyfits.open("../true_node_vals.fits")

true_node_vals = f[0].data
f.close()

try:
    f = pyfits.open("true_scene.fits")
except:
    f = pyfits.open("../true_scene.fits")

true_scene = f[0].data
f.close()

true_scene_reshape = reshape(true_scene, len(j))

print j.shape, true_scene.shape

npar = len(true_node_vals)

j = concatenate((j, diag(regularization*ones(npar, dtype=float64))))
true_scene_reshape = concatenate((true_scene_reshape, zeros(npar, dtype=float64)))

print j.shape

solved = linalg.lstsq(j, true_scene_reshape)[0]

est_model = dot(j[:-npar], solved)

save_img(true_scene - reshape(est_model, true_scene.shape), "resid.fits")

print abs(solved - true_node_vals).max()
save_img(est_model, "est_model.fits")
save_img(solved, "est_params.fits")
