import commands
from numpy import *

def randsym(sc):
    return random.random()*sc - sc/2.


def writefl(wd):
    f = open(wd + "/paramfile.txt", 'w')
    f.write("""# Arcseconds:
pixel_scale		0.05
# Arcseconds:
slice_scale		0.15
# Arcseconds:
basis_step		0.05
# Relative to FoV:
max_spline_rad		0.6

# In Angstroms:
wave		        """ + str(wv*10000) + """

n_slice			20
# Number of wavelengths to simulate. Must be large enough for n_subwave
n_wave			10
# Number of sub-sampled wavelengths:
n_subwave		61
oversample		10

# Arcseconds, not pixels!
slice_width		3

# Arcseconds, not pixels!
# Parallel to slit
xdithers		""" + str(xdithers) + """
# Perpendicular to slit
ydithers		""" + str(ydithers) + """
# Radians:
thetadithers		""" + str(thetadithers) + """
scaledithers		""" + str(scaledithers) + """

# Number of cores (max = n_slice)
processes		10

max_noise_to_signal	""" + str(0.05*include_noise) + """

galaxy_model            "load"
jacobian_instead	1
""")

raw_input("Removing old results. Proceed?")
    
xdithers = [randsym(0.05) for i in range(5)]
ydithers = [randsym(0.05) for i in range(5)]
thetadithers = [0., pi/4., pi/2., pi*3./4., pi/8.]
scaledithers = [1]*4 + [0.01]

fw = open("run_wrap.sh", 'w')

pwd = commands.getoutput("pwd")

for j in range(4):
    for include_noise in [0, 1]:
        for wv in [0.8, 0.5, 1.0]:
            wd = "test_run_ndith=%02i_nrot=%02i_wv=%.2f_%02i_noise=%i" % (len(xdithers), len(unique(thetadithers)), wv, j, include_noise)
            commands.getoutput("rm -fr " + wd)
            commands.getoutput("mkdir " + wd)
            writefl(wd)
            
            
            fw.write("cd " + pwd + "/" + wd + '\n')
            fw.write("python ../model.py paramfile.txt > log1.txt\n")
            fw.write("python ../solve.py 0 > log2.txt &\n")
fw.close()



