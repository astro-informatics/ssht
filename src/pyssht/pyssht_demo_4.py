from __future__ import print_function, unicode_literals
import numpy as np
import pyssht as ssht
import matplotlib.pyplot as plt
from pylab import cm

# % pyssht_demo_4 - Run demo4
# %
# % Demo to compute inverse and forward transform of spin function, using
# % polar interface.
# %
# % Author: Christopher G R Wallis & Jason McEwen (www.christophergrwallis.org & www.jasonmcewen.org)
# %
# % pyssht python package to perform spin spherical harmonic transforms

# Define parameters.
L = 64
Spin = 0
method = "MW_pole"

# Generate random flms (of complex signal).
flm = np.random.randn(L * L) + 1j * np.random.randn(L * L)

# Zero harmonic coefficients with el<|spin|.
ind_min = np.abs(Spin) ** 2
flm[0:ind_min] = 0.0 + 1j * 0.0

# Compute inverse then forward transform.
f, f_sp, phi_sp = ssht.inverse(flm, L, Spin=Spin, Method="MW_pole")

flm_syn = ssht.forward((f, f_sp, phi_sp), L, Spin=Spin, Method="MW_pole")

# Compute max error in harmonic space.
maxerr = np.abs(flm_syn - flm).max()
print("Max error:", maxerr)
