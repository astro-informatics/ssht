from __future__ import print_function, unicode_literals
import numpy as np
import pyssht as ssht
import matplotlib.pyplot as plt
from pylab import cm

# % pyssht_demo_3 - Run demo3
# %
# % Demo to compute inverse and forward transform of spin function, using
# % standard interface with various options.
# %
# % Author: Christopher G R Wallis & Jason McEwen (www.christophergrwallis.org & www.jasonmcewen.org)
# %
# % pyssht python package to perform spin spherical harmonic transforms

# Define parameters.
L = 64
spin = 4
methods = ["MW", "MWSS", "GL", "DH"]

# Generate random flms (of complex signal).
flm = np.random.randn(L * L) + 1j * np.random.randn(L * L)

# Zero harmonic coefficients with el<|spin|.
ind_min = np.abs(spin) ** 2
flm[0:ind_min] = 0.0 + 1j * 0.0

# Compute inverse then forward transform.
for method in methods:
    f = ssht.inverse(flm, L, Method=method, Spin=spin, Reality=False)
    flm_syn = ssht.forward(f, L, Method=method, Spin=spin, Reality=False)

    # Compute max error in harmonic space.
    maxerr = np.abs(flm_syn - flm).max()
    print("Method:", method, "\nMax error:", maxerr)
