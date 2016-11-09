import numpy as np
import pyssht as ssht
import matplotlib.pyplot as plt
from pylab import cm
import scipy.io as sio

# % pyssht_demo_5 - Run demo5
# %
# % Smooth Earth topography map by applying a Gaussian filter in harmonic
# % space.  
# %
# % The official Earth Gravitational Model EGM2008 has been publicly
# % released by the U.S. National Geospatial-Intelligence Agency (NGA) (these
# % data were downloaded from Frederik Simons' web page:
# % http://geoweb.princeton.edu/people/simons).
# %
# % Author: Christopher G R Wallis & Jason McEwen (www.christophergrwallis.org & www.jasonmcewen.org)
# %
# % pyssht python package to perform spin spherical harmonic transforms

# Define parameters.
L = 128
sigma = np.pi/L;

#% Load harmonic coefficients of Earth.
mat_contents = sio.loadmat('src/matlab/data/EGM2008_Topography_flms_L0128')

flm     = np.ascontiguousarray(mat_contents['flm'][:,0])

#% Smooth harmonic coefficients.
flm_smooth = ssht.guassian_smoothing(flm, L, sigma)

# Compute real space version of Earth.
f = ssht.inverse(flm, L, Reality=True)
f_smooth = ssht.inverse(flm_smooth, L, Reality=True)


# Plot 
f_plot, mask_array = ssht.mollweide_projection(f, L, resolution=200, rot=[0.0,np.pi,np.pi])
plt.figure()
plt.subplot(1,2,1)
imgplot = plt.imshow(f_plot,interpolation='nearest')
plt.colorbar(imgplot,fraction=0.025, pad=0.04)
plt.imshow(mask_array, interpolation='nearest', cmap=cm.gray, vmin=-1., vmax=1.)
plt.gca().set_aspect("equal")
plt.title("f")
plt.axis('off')

f_plot, mask_array = ssht.mollweide_projection(f_smooth, L, resolution=200, rot=[0.0,np.pi,np.pi])
plt.subplot(1,2,2)
imgplot = plt.imshow(f_plot,interpolation='nearest')
plt.colorbar(imgplot,fraction=0.025, pad=0.04)
plt.imshow(mask_array, interpolation='nearest', cmap=cm.gray, vmin=-1., vmax=1.)
plt.gca().set_aspect("equal")
plt.title("f smooth")
plt.axis('off')

plt.show()