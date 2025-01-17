from numpy import *
import commands
import pyfits
from scipy.interpolate import RectBivariateSpline, interp1d
import time
import multiprocessing as mp
from matplotlib import use
use("PDF")
import matplotlib.pyplot as plt
import sys
from scipy.optimize import minimize


def save_img(dat, imname, waves = None):
    commands.getoutput("rm -f " + imname)
    fitsobj = pyfits.HDUList()
    hdu = pyfits.PrimaryHDU()
    hdu.data = dat
    fitsobj.append(hdu)
    fitsobj.writeto(imname)
    fitsobj.close()

def get_params():
    f = open(sys.argv[1])
    lines = f.read().split('\n')
    f.close()

    params = {}
    for line in lines:
        if line.count("#") == 0:
            parsed = line.split(None)
            if len(parsed) > 1:
                params[parsed[0]] = eval(" ".join(parsed[1:]))
    return params
    

def load_basis_fn():
    basisfns = {}
    for item in unique(thetadithers):
        f = pyfits.open(basis_fl % item)
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

    for wave in range(n_subwave):
        the_range = arange(oversample*slice_scale/pixel_scale)
        output_pixel = array(floor((wave + the_range) / float(oversample)), dtype=int32)

        
        the_ranges = []
        for the_out in sort(unique(output_pixel)):
            inds = where(output_pixel == the_out)[0]
            
            the_ranges.append([inds[0], inds[-1]+1, the_out])
        print wave, the_ranges
        wave_to_inds.append(the_ranges)
    return wave_to_inds


def single_param_jacobian(args):
    i, basis_fns, nudgex, nudgey = args

    tmpout = zeros([len(xdithers), n_slice, slice_width/pixel_scale, n_wave], dtype=float64)
    
    for m in range(len(xdithers)):
        nodex = node_xy[i][1]*cos(thetadithers[m]) + node_xy[i][2]*sin(thetadithers[m])
        nodey = -node_xy[i][1]*sin(thetadithers[m]) + node_xy[i][2]*cos(thetadithers[m])
        
        basis_eval = basis_fns[thetadithers[m]](xvals - nodex - xdithers[m] + 0.0025 - nudgex, yvals - nodey - ydithers[m] + 0.005 - nudgey)
        
        
        for k in range(n_slice):
            for n in range(n_subwave): # Subsampled wavelength                                                                                
                unbinned_term = basis_eval[:,k*30:(k+1)*30]
                
                for item in wave_to_inds[n]:
                    tmpout[m, k, :,item[2]] += sum(unbinned_term[:,item[0]:item[1]], axis = 1)*wavelength_splines[node_xy[i][3],n]*scaledithers[m]

    return reshape(tmpout, len(xdithers)*n_slice*(slice_width/pixel_scale)*n_wave)



def make_jacobian(nudgex, nudgey):
    jacobian = pool.map(single_param_jacobian, [(i, basis_fns, nudgex, nudgey) for i in range(len(node_xy))])
    jacobian = array(jacobian)

    print "jacobian", jacobian.shape
    if len(jacobian) < 20000:
        jacobian = transpose(jacobian)
    print "jacobian", jacobian.shape

    return jacobian


def make_dither(inputs):
    k, node_vals, resid, node_xy, basis_fns, wavelength_splines = inputs
    outputs = zeros([len(xdithers), n_slice, slice_width/pixel_scale, n_wave], dtype=float64)
    gradient = zeros(len(node_xy), dtype=float64)

    if k % 10 == 0:
        print "Slice ", k
    for m in range(len(xdithers)):

        for i in range(len(node_xy)):
            nodex = node_xy[i][1]*cos(thetadithers[m]) - node_xy[i][2]*sin(thetadithers[m])
            nodey = node_xy[i][1]*sin(thetadithers[m]) + node_xy[i][2]*cos(thetadithers[m])

            basis_eval = basis_fns[thetadithers[m]](xvals - nodex - xdithers[m], yvals - nodey - ydithers[m])
            
            for n in range(n_subwave): # Subsampled wavelength
                
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
                                           
    
def get_scene(node_vals, xdithers, ydithers, thetadithers, node_xy, basis_fns, wavelength_splines, resid = None):
    #for k in range(n_slice):
    #    make_dither([k, node_vals, resid, node_xy, basis_fns])
    
    results = pool.map(make_dither, [(k, node_vals, resid, node_xy, basis_fns, wavelength_splines) for k in range(n_slice)])
    
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

        wavelength_splines.append(ifn(arange(n_subwave, dtype=float64)/(basis_step/pixel_scale * oversample))
                                  )
        plt.plot(wavelength_splines[-1])
    plt.savefig("wavelength_splines.pdf")
    plt.close()

    return node_xy, array(wavelength_splines)
    
def chi2fn(x):
    scene, NA = get_scene(x, xdithers, ydithers, thetadithers, node_xy, basis_fns, wavelength_splines)
    resid = true_scene - scene
    save_img(resid, "resid.fits")
    
    scene2, gradient = get_scene(x, xdithers, ydithers, thetadithers, node_xy, basis_fns, wavelength_splines, resid = resid)

    print sum(resid**2.), log10(sum(resid**2.))
    return sum(resid**2.), gradient



"""
Coordinates:
    Scene: Arcseconds, -1.5 to +1.5; parallel to slice first

"""

params = get_params()

pixel_scale = params["pixel_scale"]
slice_scale = params["slice_scale"]
n_slice = params["n_slice"]
slice_width = params["slice_width"]
n_wave = params["n_wave"]
basis_step = params["basis_step"]
oversample = params["oversample"]
n_subwave = params["n_subwave"]

xdithers = params["xdithers"]
ydithers = params["ydithers"]
thetadithers = params["thetadithers"]
scaledithers = params["scaledithers"]

wavelen = params["wave"]/10000. # Convert to microns

########################################## Done with setup ##########################################

assert isclose(n_subwave % (basis_step/pixel_scale * oversample), 1.0)
n_wave_spline = int(around(n_subwave/(basis_step/pixel_scale * oversample))) + 1
print "n_wave_spline", n_wave_spline

wave_to_inds = get_wave_to_inds()

node_xy, wavelength_splines = get_spline_nodes(scale = params["max_spline_rad"])
print "number of nodes", len(node_xy)

xvals = arange(-1.475, 1.5, 0.05)
yvals = arange(-1.5, 1.5, 0.005)

pool = mp.Pool(processes = params["processes"])


########################################## Done with derived quantities ##########################################


#true_node_vals = arange(len(node_xy)*n_wave_spline, dtype=float64)

if params["galaxy_model"] == "random":
    basis_fl = "../../psfs/basis_5_mas_rot=%s_spl=0.050_wv=%.2f.fits" % ("%.3f", wavelen)
    basis_fns = load_basis_fn()


    true_node_vals = random.normal(size = len(node_xy))
    save_img(true_node_vals, "true_node_vals.fits")

    t = time.time()
    true_scene, NA = get_scene(true_node_vals, xdithers, ydithers, thetadithers, node_xy, basis_fns, wavelength_splines)
    t2 = time.time()
    
    print "True scene in ", t2 - t

if params["galaxy_model"] == "points":
    basis_fl = "../../psfs/basis_5_mas_rot=%s_spl=0.050_wv=%.2f.fits" % ("%.3f", wavelen)
    basis_fns = load_basis_fn()


    true_node_vals = zeros(len(node_xy), dtype=float64)
    true_node_vals[0] = 1
    true_node_vals[-1] = 0.5

    save_img(true_node_vals, "true_node_vals.fits")

    t = time.time()
    true_scene, NA = get_scene(true_node_vals, xdithers, ydithers, thetadithers, node_xy, basis_fns, wavelength_splines)
    t2 = time.time()
    
    print "True scene in ", t2 - t

if params["galaxy_model"] == "load":
    basis_fl = "../../psfs/galaxybasis_5_mas_rot=%s_wv=%.2f.fits" % ("%.3f", wavelen)
    basis_fns = load_basis_fn()
    
    spec_fl = "../../psfs/spectrum_subsamp=10x_wv=%.2f.txt" % wavelen
    spec_oversample = loadtxt(spec_fl)
    

    n_wave_spline_tmp = n_wave_spline
    n_wave_spline = n_subwave
    wavelength_splines_tmp = wavelength_splines
    wavelength_splines = identity(n_subwave, dtype=float64)
    wave_to_inds = get_wave_to_inds()

    
    tmp_node_xy = [(0, 0, 0, i) for i in range(n_wave_spline)]
    true_node_vals = spec_oversample[:n_subwave]
    # Now, apodize
    for i in range(12):
        # i ranges from 0 to 11
        sci = i/11.
        true_node_vals[i] *= 3*sci**2. - 2.*sci**3.
        true_node_vals[-1 -i] *= 3*sci**2. - 2.*sci**3.

    print "true_node_vals", true_node_vals

    save_img(true_node_vals, "true_node_vals.fits")

    print "Here"
    t = time.time()
    
    true_scene, NA = get_scene(true_node_vals, xdithers, ydithers, thetadithers, tmp_node_xy, basis_fns, wavelength_splines)
    t2 = time.time()

    print "True scene in ", t2 - t
    print "Made scene"

    n_wave_spline = n_wave_spline_tmp
    wavelength_splines = wavelength_splines_tmp
    wave_to_inds = get_wave_to_inds()

    basis_fl = "../../psfs/basis_5_mas_rot=%s_spl=0.050_wv=%.2f.fits" % ("%.3f", wavelen)
    basis_fns = load_basis_fn()



if params["max_noise_to_signal"] > 0:
    save_img(true_scene, "noise_free_scene.fits")
    noise_scale = abs(true_scene).max()*params["max_noise_to_signal"]

    for i in range(len(scaledithers)):
        true_scene[i] += random.normal(size = true_scene[i].shape)*noise_scale*scaledithers[i]

#true_scene /= true_scene.max()
#true_scene += random.normal(size = true_scene.shape)*0.1

save_img(true_scene, "true_scene.fits")




if params["jacobian_instead"] == 1:
    print "Making original jacobian..."
    jacobian = make_jacobian(nudgex = 0, nudgey = 0)
    save_img(jacobian, "jacobian0.fits")

    true_scene_reshape = reshape(true_scene, len(jacobian))

    print jacobian.shape, true_scene.shape

    npar = len(jacobian[0])

    print "Solving original jacobian..."

    solved = linalg.lstsq(jacobian, true_scene_reshape)[0]

    est_model0 = dot(jacobian, solved)
    est_model = reshape(est_model0, true_scene.shape)

    resid = true_scene - est_model
    save_img(resid, "resid0.fits")

    print "Making dx jacobian..."
    jacobiandx = make_jacobian(nudgex = 0.01, nudgey = 0)
    print "Making dy jacobian..."
    jacobiandy = make_jacobian(nudgex = 0.0, nudgey = 0.01)

    est_modeldx = dot(jacobiandx, solved)
    est_modeldy = dot(jacobiandy, solved)
    
    jdx = (est_modeldx - est_model0)/0.01
    jdy = (est_modeldy - est_model0)/0.01
    
    jcomb = zeros([len(est_model0), 2], dtype=float64)
    jcomb[:,0] = jdx
    jcomb[:,1] = jdy

    print "Solving for dx and dy..."

    solvedcomb = linalg.lstsq(jcomb, true_scene_reshape - est_model0)[0]
    print "Solved dx, dy ", solvedcomb
    
    print "Making final jacobian..."

    jacobian = make_jacobian(nudgex = solvedcomb[0], nudgey = solvedcomb[1])
    save_img(jacobian, "jacobian.fits")

    print "Doing final solve..."

    solved = linalg.lstsq(jacobian, true_scene_reshape)[0]
    est_model0 = dot(jacobian, solved)
    est_model = reshape(est_model0, true_scene.shape)
    resid = true_scene - est_model
    save_img(resid, "resid.fits")




    #true_scene2 = dot(jacobian, true_node_vals)
    #true_scene2 = reshape(true_scene2, true_scene.shape)
    #print abs(true_scene - true_scene2).max()
    stop_now




opt_res = minimize(chi2fn, x0 = zeros(len(node_xy)), method="L-BFGS-B", jac = True, options = dict(disp = True, maxcor = 30))
#opt_res = minimize(chi2fn, x0 = zeros(len(node_xy)), method="BFGS", jac = True, options = dict(disp = True))
#opt_res = minimize(chi2fn, x0 = zeros(len(node_xy)), method="CG", jac = True, options = dict(disp = True))
#opt_res = minimize(chi2fn, x0 = zeros(len(node_xy)), method="Newton-CG", jac = True, options = dict(disp = True))
