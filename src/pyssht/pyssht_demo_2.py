from __future__ import print_function
import numpy as np
import pyssht as ssht
import matplotlib.pyplot as plt
from pylab import cm

# % pyssht_demo_2 - Run demo2
# %
# % Simple demo to compute inverse and forward transform of real scalar
# % function, using simplest interface with default options.
# %
# % Author: Christopher G R Wallis & Jason McEwen (www.christophergrwallis.org & www.jasonmcewen.org)
# %
# % pyssht python package to perform spin spherical harmonic transforms

# Define parameters.
L = 64

# Generate random flms (of real signal).
flm = np.zeros((L * L), dtype=complex)

# Impose reality on flms.
for el in range(L):
    m = 0
    ind = ssht.elm2ind(el, m)
    flm[ind] = np.random.randn()
    for m in range(1, el + 1):
        ind_pm = ssht.elm2ind(el, m)
        ind_nm = ssht.elm2ind(el, -m)
        flm[ind_pm] = np.random.randn() + 1j * np.random.randn()
        flm[ind_nm] = (-1) ** m * np.conj(flm[ind_pm])

# Compute inverse then forward transform.
f = ssht.inverse(flm, L, Reality=True)
flm_syn = ssht.forward(f, L, Reality=True)

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
