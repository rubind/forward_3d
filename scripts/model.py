from numpy import *
import commands
import pyfits
from scipy.interpolate import RectBivariateSpline, interp1d
import time
import multiprocessing as mp
from matplotlib import use
use("PDF")
import matplotlib.pyplot as plt


def save_img(dat, imname, waves = None):
    commands.getoutput("rm -f " + imname)
    fitsobj = pyfits.HDUList()
    hdu = pyfits.PrimaryHDU()
    hdu.data = dat
    fitsobj.append(hdu)
    fitsobj.writeto(imname)
    fitsobj.close()

def load_basis_fn():
    basisfns = {}
    for item in unique(thetadithers):
        f = pyfits.open("../psfs/basis_5_mas_rot=%.3f.fits" % item)
        basis = f[0].data
        f.close()
        vals0 = arange(len(basis), dtype=float64)
        vals1 = arange(len(basis), dtype=float64)
        
        inds = where(basis == basis.max())
        vals0 -= vals0[inds[0]]
        vals1 -= vals1[inds[0]]
        
        vals0 *= pixel_scale/float(oversample)
        vals1 *= pixel_scale/float(oversample)
        
        basisfns[item] = RectBivariateSpline(vals0, vals1, basis, kx = 3, ky = 3, s = 0)

    return basisfns

    
def bin_to_pixels(to_bin, the_output, wave, scale = 1.):
    assert to_bin.shape == (slice_width/pixel_scale, oversample*slice_scale/pixel_scale) # E.g., 60 by 30

    wave_frac = wave % oversample
    wave_int = wave - wave_frac

    to_bin_pad = pad(to_bin, [(0,0), (wave_frac, (oversample - wave_frac) % oversample)], 'constant', constant_values = 0)

    padded_size = to_bin_pad.shape[1]
    
    binned_output = 0.
    for i in range(oversample):
        binned_output += to_bin_pad[:,i::oversample]

    the_output[:, wave_int/oversample: (wave_int + padded_size)/oversample] += binned_output*scale

    return the_output

def get_wave_to_inds():
    wave_to_inds = []

    for wave in range(hard_code):
        the_range = arange(oversample*slice_scale/pixel_scale)
        output_pixel = array(floor((wave + the_range) / float(oversample)), dtype=int32)

        
        the_ranges = []
        for the_out in sort(unique(output_pixel)):
            inds = where(output_pixel == the_out)[0]
            
            the_ranges.append([inds[0], inds[-1]+1, the_out])
        print wave, the_ranges
        wave_to_inds.append(the_ranges)
    return wave_to_inds


def make_dither(inputs):
    k, node_vals, resid = inputs
    outputs = zeros([len(xdithers), n_slice, slice_width/pixel_scale, n_wave], dtype=float64)
    gradient = zeros(len(node_xy), dtype=float64)

    if k % 10 == 0:
        print "Slice ", k
    for m in range(len(xdithers)):

        for i in range(len(node_xy)):
            nodex = node_xy[i][1]*cos(thetadithers[m]) - node_xy[i][2]*sin(thetadithers[m])
            nodey = node_xy[i][1]*sin(thetadithers[m]) + node_xy[i][2]*cos(thetadithers[m])

            basis_eval = basis_fns[thetadithers[m]](xvals - nodex - xdithers[m], yvals - nodey - ydithers[m])
            
            for n in range(hard_code): # Subsampled wavelength
                
                #this_term = zeros([slice_width/pixel_scale, n_wave])
                #this_term = bin_to_pixels(basis_eval[:,k*30:(k+1)*30], this_term, wave = n, scale = 1.0)
                unbinned_term = basis_eval[:,k*30:(k+1)*30]
                whole_term = outputs[m,k]*0.
                
                for item in wave_to_inds[n]:
                    whole_term[:,item[2]] = sum(unbinned_term[:,item[0]:item[1]], axis = 1)

                outputs[m, k] += whole_term*node_vals[i]*wavelength_splines[node_xy[i][3],n]*scaledithers[m]
                if resid != None:
                    gradient[i] += -2.*sum(whole_term*resid[m,k]*wavelength_splines[node_xy[i][3],n]*scaledithers[m])

    return outputs, gradient
                                           
    
def get_scene(node_vals, xdithers, ydithers, thetadithers, resid = None):
    results = pool.map(make_dither, [(k, node_vals, resid) for k in range(n_slice)])
        
    outputs = 0.
    gradient = 0.

    for k in range(n_slice):
        outputs += results[k][0]
        if resid != None:
            gradient += results[k][1]

    return outputs, gradient

def get_spline_nodes(scale = 0.9):
    node_xy = []
    
    for dx in arange(-basis_step*50., basis_step*50., basis_step):
        for dy in arange(-basis_step*50., basis_step*50., basis_step):
            r = sqrt(dx**2. + dy**2.)
            if r < slice_width*scale/2.:
                for i in range(n_wave_spline):
                    node_xy.append((r, dx, dy, i))

    node_xy.sort()

    wavelength_splines = []
    for i in range(n_wave_spline):
        xvals = arange(n_wave_spline, dtype=float64)
        yvals = xvals*0.
        
        yvals[i] = 1.

        ifn = interp1d(xvals, yvals, kind = 'cubic')

        wavelength_splines.append(ifn(arange(hard_code, dtype=float64)/(basis_step/pixel_scale * oversample))
                                  )
        plt.plot(wavelength_splines[-1])
    plt.savefig("wavelength_splines.pdf")
    plt.close()

    return node_xy, array(wavelength_splines)
    



"""
Coordinates:
    Scene: Arcseconds, -1.5 to +1.5; parallel to slice first

"""

pixel_scale = 0.05 # Arcseconds
slice_scale = 0.15 # Arcseconds
n_slice = 20
slice_width = 3 # Arcseconds, not pixels!
n_wave = 7
basis_step = 0.035 # Arcseconds
oversample = 10
hard_code = 36 # Number of sub-sampled wavelengths

xdithers = [0., 0.025, 0, 0.025]*2 + [0.02]
ydithers = [0., 0., 0.025, 0.025]*2 + [0.02]
thetadithers = [0]*4 + [pi/2.]*4 + [pi/4.]
scaledithers = [1.]*8 + [0.01]

########################################## Done with setup ##########################################

assert isclose(hard_code % (basis_step/pixel_scale * oversample), 1.0)
n_wave_spline = int(around(hard_code/(basis_step/pixel_scale * oversample))) + 1
print "n_wave_spline", n_wave_spline

wave_to_inds = get_wave_to_inds()
basis_fns = load_basis_fn()
node_xy, wavelength_splines = get_spline_nodes(scale = 0.5)

xvals = arange(-1.475, 1.5, 0.05)
yvals = arange(-1.5, 1.5, 0.005)

pool = mp.Pool(processes = 10)

########################################## Done with derived quantities ##########################################


#true_node_vals = arange(len(node_xy)*n_wave_spline, dtype=float64)

true_node_vals = random.normal(size = len(node_xy))



t = time.time()
true_scene, NA = get_scene(true_node_vals, xdithers, ydithers, thetadithers)
t2 = time.time()

print "True scene in ", t2 - t

#true_scene /= true_scene.max()
#true_scene += random.normal(size = true_scene.shape)*0.1

save_img(true_scene, "true_scene.fits")


def chi2fn(x):
    scene, NA = get_scene(x, xdithers, ydithers, thetadithers)
    resid = true_scene - scene
    save_img(resid, "resid.fits")
    
    scene2, gradient = get_scene(x, xdithers, ydithers, thetadithers, resid = resid)
    return sum(resid**2.), gradient

from scipy.optimize import minimize

opt_res = minimize(chi2fn, x0 = zeros(len(node_xy)), method="L-BFGS-B", jac = True, options = dict(disp = True, maxcor = 30))
