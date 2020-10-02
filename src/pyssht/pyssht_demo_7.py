import numpy as np
import pyssht as ssht
import matplotlib.pyplot as plt
from pylab import cm
import scipy.io as sio

# % pyssht_demo_7 - Run demo7
# %
# % Plot spherical harmonic function on the sphere, then
# % rotate it and plot the rotation, too.
# %
# % Author: Christopher G R Wallis & Jason McEwen (www.christophergrwallis.org & www.jasonmcewen.org)
# %
# % pyssht python package to perform spin spherical harmonic transforms

# Define parameters.
L = 64
el = 4
m = 2
gamma = np.pi / 2
beta = np.pi / 4
alpha = -np.pi / 2

# Generate spherical harmonics.
flm = np.zeros((L * L), dtype=complex)
ind = ssht.elm2ind(el, m)
flm[ind] = 1.0

# Compute function on the sphere.
f = ssht.inverse(flm, L)

# Rotate spherical harmonic
flm_rot = ssht.rotate_flms(flm, alpha, beta, gamma, L)

# Compute rotated function on the sphere.
f_rot = ssht.inverse(flm_rot, L)

# Plot
f_plot, mask_array, f_plot_imag, mask_array_imag = ssht.mollweide_projection(
    f, L, resolution=200, rot=[0.0, np.pi, np.pi]
)
plt.figure()
plt.subplot(1, 2, 1)
imgplot = plt.imshow(f_plot, interpolation="nearest")
plt.colorbar(imgplot, fraction=0.025, pad=0.04)
plt.imshow(mask_array, interpolation="nearest", cmap=cm.gray, vmin=-1.0, vmax=1.0)
plt.gca().set_aspect("equal")
plt.title("f")
plt.axis("off")

f_plot, mask_array, f_plot_imag, mask_array_imag = ssht.mollweide_projection(
    f_rot, L, resolution=200, rot=[0.0, np.pi, np.pi]
)
plt.subplot(1, 2, 2)
imgplot = plt.imshow(f_plot, interpolation="nearest")
plt.colorbar(imgplot, fraction=0.025, pad=0.04)
plt.imshow(mask_array, interpolation="nearest", cmap=cm.gray, vmin=-1.0, vmax=1.0)
plt.gca().set_aspect("equal")
plt.title("f rot")
plt.axis("off")

plt.show()
