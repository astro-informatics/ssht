from __future__ import print_function, unicode_literals
import numpy as np
import pyssht as ssht
import matplotlib.pyplot as plt
from pylab import cm
import scipy.io as sio

# % pyssht_demo_6 - Run demo6
# %
# % Integrate a band-limited function on the sphere using the symmetrised
# % quadrature weights.
# %
# % Author: Christopher G R Wallis & Jason McEwen (www.christophergrwallis.org & www.jasonmcewen.org)
# %
# % pyssht python package to perform spin spherical harmonic transforms


def mw_weights(m):
    # mw_weights - Compute MW weights
    #
    # Compute the MW weights.
    #
    # Default usage is given by
    #
    #   w = mw_weights(m)
    #
    # where m is the weight index from -(L-1):(L-1) and w is the corresponding
    # weight.
    #
    # Author: Jason McEwen (www.jasonmcewen.org)

    if m == 1:
        w = 1j * np.pi / 2
    elif m == -1:
        w = -1j * np.pi / 2
    elif m % 2 == 0:
        w = 2.0 / (1.0 - m * m)
    else:
        w = 0

    return w


# Define parameters.
method = "MW"  # cant choose different ones
L = 4
reality = False
renorm_plot = True

# Generate random flms (of real signal).
flm = np.zeros((L * L), dtype=complex)

# Impose reality on flms.
if reality:
    for el in range(L):
        m = 0
        ind = ssht.elm2ind(el, m)
        flm[ind] = np.random.randn()
        for m in range(1, el + 1):
            ind_pm = ssht.elm2ind(el, m)
            ind_nm = ssht.elm2ind(el, -m)
            flm[ind_pm] = np.random.randn() + 1j * np.random.randn()
            flm[ind_nm] = (-1) ** m * np.conj(flm[ind_pm])
else:
    flm = np.random.randn(L * L) + 1j * np.random.randn(L * L)

# Compute weights.
w = np.zeros(2 * L - 1, dtype=complex)
wr = np.zeros(2 * L - 1, dtype=complex)
for i, m in enumerate(range(-(L - 1), L)):
    w[i] = mw_weights(m)
    wr[i] = mw_weights(m) * np.exp(-1j * m * np.pi / (2 * L - 1))  # apply phase shift

# Compute weights as function of theta.
wr = (np.fft.fft(np.fft.ifftshift(wr)) * 2 * np.pi / (2 * L - 1) ** 2).real

# Compute symmetrised quadrature weights defined on sphere.
q = wr[0:L]

for i, j in enumerate(range(2 * L - 2, L - 1, -1)):
    q[i] = q[i] + wr[j]

# Integral of function given by rescaled (el,m)=(0,0) harmonic coefficient.
I0 = flm[0] * np.sqrt(4 * np.pi)
print("Integration using input flm:", I0)

# Integrate function on sphere using all points.
f = ssht.inverse(flm, L, Method=method, Reality=reality)
Q = np.outer(q, np.ones(2 * L - 1))
I1 = sum(sum(Q[:] * f[:]))
print("Integration using all points:", I1)

# Integration function on sphere using L points in phi.
# (Zero-pad to evalue the band-limited function on the correct L points.)
ntheta, nphi = ssht.sample_shape(L, Method=method)

f_up = np.zeros((L, 2 * L), dtype=complex)
f_down = np.zeros((L, int((2 * L) / 2)), dtype=complex)
for r in range(ntheta):
    dum = np.fft.fftshift(np.fft.fft(f[r, :]))
    dum = np.insert(dum, 0, 0.0)
    f_up[r, :] = np.fft.ifft(np.fft.ifftshift(dum)) / (2 * L - 1) * (2 * L)

for i in range(int((nphi + 1) / 2)):
    f_down[:, i] = f_up[:, i * 2]
Q_down = np.outer(q, np.ones(L))
I2 = sum(sum(Q_down * f_down)) * (2 * L - 1) / L
print("Integration using L points in phi:", I2)

# Compute integration errors.
I1_err = abs(I1 - I0)
I2_err = abs(I2 - I0)
print("Error:", I1_err, I2_err)
