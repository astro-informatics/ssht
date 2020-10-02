from __future__ import print_function
import numpy as np
import pyssht as ssht
import matplotlib.pyplot as plt
from pylab import cm

# % pyssht_demo_1 - Run demo1
# %
# % Simple demo to compute inverse and forward transform of complex scalar
# % function, using simplest interface with default options.
# %
# % Author: Christopher G R Wallis & Jason McEwen (www.christophergrwallis.org & www.jasonmcewen.org)
# %
# % pyssht python package to perform spin spherical harmonic transforms

# Define parameters.
L = 64

# Generate random flms (of complex signal).
flm = np.random.randn(L * L) + 1j * np.random.randn(L * L)

# Compute inverse then forward transform.
f = ssht.inverse(flm, L)
flm_syn = ssht.forward(f, L)

# Compute max error in harmonic space.
maxerr = np.abs(flm_syn - flm).max()
print("Max error:", maxerr)

# Plot function on sphere using mollweide projection
f_plot, mask_array = ssht.mollweide_projection(np.abs(f), L, resolution=200)
plt.figure()
imgplot = plt.imshow(f_plot, interpolation="nearest")
plt.colorbar(imgplot, fraction=0.025, pad=0.04)
plt.imshow(mask_array, interpolation="nearest", cmap=cm.gray, vmin=-1.0, vmax=1.0)
plt.gca().set_aspect("equal")
plt.title("|f|")
plt.axis("off")
plt.show()

plt.figure()
