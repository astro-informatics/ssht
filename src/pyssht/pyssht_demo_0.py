from __future__ import unicode_literals
import numpy as np
import pyssht as ssht
import matplotlib.pyplot as plt
from pylab import cm

# % pyssht_demo_0 - Run demo0
# %
# % Plot spherical harmonic functions on the sphere using othographic projection
# %
# %
# % Author: Christopher G R Wallis & Jason McEwen (www.christophergrwallis.org & www.jasonmcewen.org)
# %
# % pyssht python package to perform spin spherical harmonic transforms

# % Define parameters.
L = 64
el = 4
m = 2

# Generate spherical harmonics.
flm = np.zeros(L * L, dtype=complex)
ind = ssht.elm2ind(el, m)
flm[ind] = 1.0

# Compute function on the sphere.
f = ssht.inverse(flm, L)

# Plot function on sphere.
(
    f_north_real,
    mask_north_real,
    f_south_real,
    mask_south_real,
    f_north_imag,
    mask_north_imag,
    f_south_imag,
    mask_south_imag,
) = ssht.polar_projection(
    f, L, resolution=200, Method="MW", Projection="OP", rot=[0.0, np.pi / 2, 0.0]
)

plt.figure(1)
plt.subplot(2, 2, 1)
plt.imshow(f_north_real, interpolation="nearest")
plt.imshow(mask_north_real, interpolation="nearest", cmap=cm.gray, vmin=-1.0, vmax=1.0)
plt.axis("off")
plt.title("f East Real")

plt.subplot(2, 2, 2)
plt.imshow(f_north_imag, interpolation="nearest")
plt.imshow(mask_north_imag, interpolation="nearest", cmap=cm.gray, vmin=-1.0, vmax=1.0)
plt.axis("off")
plt.title("f East Imaginary")

plt.figure(1)
plt.subplot(2, 2, 3)
plt.imshow(f_south_real, interpolation="nearest")
plt.imshow(mask_south_real, interpolation="nearest", cmap=cm.gray, vmin=-1.0, vmax=1.0)
plt.axis("off")
plt.title("f West Real")

plt.subplot(2, 2, 4)
plt.imshow(f_south_imag, interpolation="nearest")
plt.imshow(mask_south_imag, interpolation="nearest", cmap=cm.gray, vmin=-1.0, vmax=1.0)
plt.axis("off")
plt.title("f West Imaginary")

plt.show()
