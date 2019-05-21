from __future__ import print_function, unicode_literals
import numpy as np
import pyssht as ssht
import matplotlib.pyplot as plt
import unittest

# s2c test
L = 256
thetas, phis = ssht.sample_positions(L, Grid=True)

f = np.zeros((L, 2 * L - 1), dtype=np.float_) + np.random.randn(L, 2 * L - 1)
# ssht.plot_sphere(phis, L,Parametric=False, Output_File='test.pdf',Show=False, Color_Bar=True, Units='Radians')

(x, y, z) = ssht.s2_to_cart(thetas, phis)
(x, y, z) = ssht.spherical_to_cart(np.ones(thetas.shape), thetas, phis)

# test rotations
flm = ssht.forward(phis, L, Reality=True)
f = ssht.inverse(flm, L, Reality=True)

flm_prime = ssht.rotate_flms(flm, np.pi / 4, np.pi / 4, np.pi / 4, L)
f_prime = ssht.inverse(flm_prime, L, Reality=True)

# ssht.plot_sphere(f, L,Parametric=False, Output_File='test_phi_sphere.pdf',Show=False, Color_Bar=True, Units='Radians')
# ssht.plot_sphere(f_prime, L,Parametric=False, Output_File='test_phi_rot_sphere.pdf',Show=False, Color_Bar=True, Units='Radians')

# plot = ssht.plot_mollweide(f, L, Close=True)
# plt.show()
# plot2 = ssht.plot_mollweide(f_prime, L, Close=True)
# plt.show()

# transform tests
L = 128
Spin = 0

flm = np.random.randn(L * L) + 1j * np.random.randn(L * L)

f = ssht.inverse(flm, L)

flm_rec = ssht.forward(f, L)

f_rec = ssht.inverse(flm_rec, L)

if (np.mean(np.abs(flm_rec - flm)) < 1E-14):
    print('MW complex transform test passed, max error:',
          np.max(np.abs(flm_rec - flm)))
else:
    print('MW complex transform ***NOT*** test passed, max error:',
          np.max(np.abs(flm_rec - flm)))

f_prime = np.zeros((L, 2 * L - 1), dtype=np.complex_) + \
    np.random.randn(L, 2 * L - 1) + 1j * np.random.randn(L, 2 * L - 1)
flm_prime = ssht.inverse_adjoint(f_prime, L)

error = np.abs(np.vdot(flm_prime, flm) - np.vdot(f_prime, f))

flm_prime = np.random.randn(L * L) + 1j * np.random.randn(L * L)
f_prime = ssht.forward_adjoint(flm_prime, L)

error += np.abs(np.vdot(flm_prime, flm) - np.vdot(f_prime, f))

if (error < 1E-5):
    print('         and the adjoints, error:', error)
else:
    print('         ****the adjoints failed****', error)

L = 128
Spin = 2

f = np.zeros((L, 2 * L - 1), dtype=np.complex_) + \
    np.random.randn(L, 2 * L - 1) + 1j * np.random.randn(L, 2 * L - 1)

flm = ssht.forward(f, L, Spin=2)

f = ssht.inverse(flm, L, Spin=2)

flm = ssht.forward(f, L, Spin=2)

f_rec = ssht.inverse(flm, L, Spin=2)

if (np.mean(np.abs(f_rec - f)) < 1E-14):
    print('MW spin complex transform test passed, max error:', np.max(
        np.abs(f_rec - f)))
else:
    print('MW spin complex transform ***NOT*** test passed, max error:', np.max(
        np.abs(f_rec - f)))

f_prime = np.zeros((L, 2 * L - 1), dtype=np.complex_) + \
    np.random.randn(L, 2 * L - 1) + 1j * np.random.randn(L, 2 * L - 1)
flm_prime = ssht.inverse_adjoint(f_prime, L, Spin=Spin)

error = np.abs(np.vdot(flm_prime, flm) - np.vdot(f_prime, f))

flm_prime = np.random.randn(L * L) + 1j * np.random.randn(L * L)
flm_prime[0:Spin * Spin] = 0.0
f_prime = ssht.forward_adjoint(flm_prime, L, Spin=Spin)

error += np.abs(np.vdot(flm_prime, flm) - np.vdot(f_prime, f))

if (error < 1E-5):
    print('         and the adjoints, error:', error)
else:
    print('         ****the adjoints failed****', error)

L = 128
spin = 0

f = np.zeros((L, 2 * L - 1), dtype=np.float_) + np.random.randn(L, 2 * L - 1)

flm = ssht.forward(f, L, Reality=True)

f = ssht.inverse(flm, L, Reality=True)

flm = ssht.forward(f, L, Reality=True)

f_rec = ssht.inverse(flm, L, Reality=True)

if (np.mean(np.abs(f_rec - f)) < 1E-14):
    print('MW real transform test passed, max error:', np.max(np.abs(f_rec - f)))
else:
    print('MW real transform ***NOT*** test passed, max error:', np.max(
        np.abs(f_rec - f)))

f_prime = np.zeros((L, 2 * L - 1), dtype=np.float_) + \
    np.random.randn(L, 2 * L - 1)

flm_prime = ssht.inverse_adjoint(f_prime, L, Reality=True)

error = np.abs(np.vdot(flm_prime, flm) - np.vdot(f_prime, f))

f_prime = np.zeros((L, 2 * L - 1), dtype=np.float_) + \
    np.random.randn(L, 2 * L - 1)
flm_prime = ssht.forward(f_prime, L, Reality=True)

f_prime = ssht.forward_adjoint(flm_prime, L, Reality=True)

error += np.abs(np.vdot(flm_prime, flm) - np.vdot(f_prime, f))

if (error < 1E-5):
    print('         and the adjoints, error:', error)
else:
    print('         ****the adjoints failed****', error)

L = 128
Spin = 0

flm = np.random.randn(L * L) + 1j * np.random.randn(L * L)

f, f_sp, phi_sp = ssht.inverse(flm, L, Method='MW_pole')

flm_rec = ssht.forward((f, f_sp, phi_sp), L, Method='MW_pole')

f_rec, f_sp_rec, phi_sp_rec = ssht.inverse(flm, L, Method='MW_pole')

if (np.mean(np.abs(flm_rec - flm)) < 1E-14):
    print('MW_pole complex transform test passed, max error:', np.max(
        np.abs(flm_rec - flm)))
else:
    print('MW_pole complex transform ***NOT*** test passed, max error:', np.max(
        np.abs(flm_rec - flm)))

L = 128
Spin = 2

flm = np.random.randn(L * L) + 1j * np.random.randn(L * L)
flm[0:4] = 0.0 + 0.0j

f, f_sp, phi_sp = ssht.inverse(flm, L, Spin=Spin, Method='MW_pole')

flm_rec = ssht.forward((f, f_sp, phi_sp), L, Spin=Spin, Method='MW_pole')

f_rec, f_sp_rec, phi_sp_rec = ssht.inverse(
    flm, L, Spin - Spin, Method='MW_pole')

if (np.mean(np.abs(flm_rec - flm)) < 1E-14):
    print('MW_pole spin complex transform test passed, max error:', np.max(
        np.abs(flm_rec - flm)))
else:
    print('MW_pole spin complex transform ***NOT*** test passed, max error:',
          np.max(np.abs(flm_rec - flm)))

L = 128
Spin = 0

f = np.random.randn(L - 1, 2 * L)

f_sp = np.random.randn()

flm = ssht.forward((f, f_sp), L, Method='MW_pole', Reality=True)

f, f_sp = ssht.inverse(flm, L, Method='MW_pole', Reality=True)

flm = ssht.forward((f, f_sp), L, Method='MW_pole', Reality=True)

f_rec, f_sp_rec = ssht.inverse(flm, L, Method='MW_pole', Reality=True)

if (np.mean(np.abs(f_rec - f)) < 1E-14 and np.abs(f_sp_rec - f_sp) < 1E-12):
    print('MW_pole real transform test passed, max error:', np.max(
        np.abs(f_rec - f)))
else:
    print('MW_pole real transform ***NOT*** test passed, max error:', np.max(
        np.abs(f_rec - f)))

L = 128
Spin = 0

flm = np.random.randn(L * L) + 1j * np.random.randn(L * L)

f = ssht.inverse(flm, L, Method='MWSS')

flm_rec = ssht.forward(f, L, Method='MWSS')

f_rec = ssht.inverse(flm, L, Method='MWSS')

if (np.mean(np.abs(flm_rec - flm)) < 1E-14):
    print('MWSS complex transform test passed, max error:', np.max(
        np.abs(flm_rec - flm)))
else:
    print('MWSS complex transform ***NOT*** test passed, max error:', np.max(
        np.abs(flm_rec - flm)))

n_theta, n_phi = ssht.sample_shape(L, Method='MWSS')

f_prime = np.zeros(ssht.sample_shape(L, Method='MWSS'), dtype=np.complex_) + \
    np.random.randn(n_theta, n_phi) + 1j * np.random.randn(n_theta, n_phi)

flm_prime = ssht.inverse_adjoint(f_prime, L, Method='MWSS')

error = np.abs(np.vdot(flm_prime, flm) - np.vdot(f_prime, f))

flm_prime = np.random.randn(L * L) + 1j * np.random.randn(L * L)
f_prime = ssht.forward_adjoint(flm_prime, L, Method='MWSS')

error += np.abs(np.vdot(flm_prime, flm) - np.vdot(f_prime, f))

if (error < 1E-5):
    print('         and the adjoints, error:', error)
else:
    print('         ****the adjoints failed****', error)

L = 128
Spin = 2

flm = np.random.randn(L * L) + 1j * np.random.randn(L * L)
flm[0:4] = 0.0 + 0.0j

f = ssht.inverse(flm, L, Spin=Spin, Method='MWSS')

flm_rec = ssht.forward(f, L, Spin=Spin, Method='MWSS')

f_rec = ssht.inverse(flm, L, Spin=Spin, Method='MWSS')

if (np.mean(np.abs(flm_rec - flm)) < 1E-14 and np.abs(f_sp_rec - f_sp) < 1E-12 and np.abs(phi_sp_rec - phi_sp) < 1E-12):
    print('MWSS spin complex transform test passed, max error:', np.max(
        np.abs(flm_rec - flm)))
else:
    print('MWSS spin complex transform ***NOT*** test passed, max error:',
          np.max(np.abs(flm_rec - flm)))

f_prime = np.zeros((n_theta, n_phi), dtype=np.complex_) + \
    np.random.randn(n_theta, n_phi) + 1j * np.random.randn(n_theta, n_phi)
flm_prime = ssht.inverse_adjoint(f_prime, L, Spin=Spin, Method='MWSS')

error = np.abs(np.vdot(flm_prime, flm) - np.vdot(f_prime, f))

flm_prime = np.random.randn(L * L) + 1j * np.random.randn(L * L)
flm_prime[0:Spin * Spin] = 0.0
f_prime = ssht.forward_adjoint(flm_prime, L, Spin=Spin, Method='MWSS')

error += np.abs(np.vdot(flm_prime, flm) - np.vdot(f_prime, f))

if (error < 1E-5):
    print('         and the adjoints, error:', error)
else:
    print('         ****the adjoints failed****', error)

L = 128
Spin = 0

f = np.random.randn(L + 1, 2 * L)

flm = ssht.forward(f, L, Method='MWSS', Reality=True)

f = ssht.inverse(flm, L, Method='MWSS', Reality=True)

flm_rec = ssht.forward(f, L, Method='MWSS', Reality=True)

f_rec = ssht.inverse(flm_rec, L, Method='MWSS', Reality=True)

if (np.mean(np.abs(f_rec - f)) < 1E-14):
    print('MWSS real transform test passed, max error:', np.max(
        np.abs(f_rec - f)))
else:
    print('MWSS real transform ***NOT*** test passed, max error:', np.max(
        np.abs(f_rec - f)))

f_prime = np.zeros((n_theta, n_phi), dtype=np.float_) + \
    np.random.randn(n_theta, n_phi)

flm_prime = ssht.inverse_adjoint(f_prime, L, Reality=True, Method='MWSS')

error = np.abs(np.vdot(flm_prime, flm) - np.vdot(f_prime, f))

f_prime = np.random.randn(n_theta, n_phi)
flm_prime = ssht.forward(f_prime, L, Reality=True)

f_prime = ssht.forward_adjoint(flm_prime, L, Reality=True, Method='MWSS')

error += np.abs(np.vdot(flm_prime, flm) - np.vdot(f_prime, f))

if (error < 1E-5):
    print('         and the adjoints, error:', error)
else:
    print('         ****the adjoints failed****', error)

L = 128
Spin = 0

f = np.zeros((2 * L, 2 * L - 1), dtype=np.complex_) + \
    np.random.randn(2 * L, 2 * L - 1) + 1j * np.random.randn(2 * L, 2 * L - 1)

flm = ssht.forward(f, L, Method='DH')

f = ssht.inverse(flm, L, Method='DH')

flm = ssht.forward(f, L, Method='DH')

f_rec = ssht.inverse(flm, L, Method='DH')

if (np.mean(np.abs(f_rec - f)) < 1E-10):
    print('DH complex transform test passed, max error:', np.max(
        np.abs(f_rec - f)))
else:
    print('DH complex transform ***NOT*** test passed, max error:', np.max(
        np.abs(f_rec - f)))

L = 128
Spin = 2

f = np.zeros((2 * L, 2 * L - 1), dtype=np.complex_) + \
    np.random.randn(2 * L, 2 * L - 1) + 1j * np.random.randn(2 * L, 2 * L - 1)

flm = ssht.forward(f, L, Spin=Spin, Method='DH')

f = ssht.inverse(flm, L, Spin=Spin, Method='DH')

flm = ssht.forward(f, L, Spin=Spin, Method='DH')

f_rec = ssht.inverse(flm, L, Spin=Spin, Method='DH')

if (np.mean(np.abs(f_rec - f)) < 1E-10):
    print('DH spin complex transform test passed, max error:', np.max(
        np.abs(f_rec - f)))
else:
    print('DH spin complex transform ***NOT*** test passed, max error:',
          np.max(np.abs(f_rec - f)))

L = 128
spin = 0

f = np.zeros((2 * L, 2 * L - 1), dtype=np.float_) + \
    np.random.randn(2 * L, 2 * L - 1)

flm = ssht.forward(f, L, Method='DH', Reality=True)

f = ssht.inverse(flm, L, Method='DH', Reality=True)

flm = ssht.forward(f, L, Method='DH', Reality=True)

f_rec = ssht.inverse(flm, L, Method='DH', Reality=True)

if (np.mean(np.abs(f_rec - f)) < 1E-10):
    print('DH real transform test passed, max error:', np.max(np.abs(f_rec - f)))
else:
    print('DH real transform ***NOT*** test passed, max error:', np.max(
        np.abs(f_rec - f)))

L = 128
Spin = 0

f = np.zeros((L, 2 * L - 1), dtype=np.complex_) + \
    np.random.randn(L, 2 * L - 1) + 1j * np.random.randn(L, 2 * L - 1)

flm = ssht.forward(f, L, Method='GL')

f = ssht.inverse(flm, L, Method='GL')

flm = ssht.forward(f, L, Method='GL')

f_rec = ssht.inverse(flm, L, Method='GL')

if (np.mean(np.abs(f_rec - f)) < 1E-10):
    print('GL complex transform test passed, max error:', np.max(
        np.abs(f_rec - f)))
else:
    print('GL complex transform ***NOT*** test passed, max error:', np.max(
        np.abs(f_rec - f)))

L = 128
Spin = 2

f = np.zeros((L, 2 * L - 1), dtype=np.complex_) + \
    np.random.randn(L, 2 * L - 1) + 1j * np.random.randn(L, 2 * L - 1)

flm = ssht.forward(f, L, Spin=Spin, Method='GL')

f = ssht.inverse(flm, L, Spin=Spin, Method='GL')

flm = ssht.forward(f, L, Spin=Spin, Method='GL')

f_rec = ssht.inverse(flm, L, Spin=Spin, Method='GL')

if (np.mean(np.abs(f_rec - f)) < 1E-10):
    print('GL spin complex transform test passed, max error:', np.max(
        np.abs(f_rec - f)))
else:
    print('GL spin complex transform ***NOT*** test passed, max error:',
          np.max(np.abs(f_rec - f)))

L = 128
spin = 0

f = np.zeros((L, 2 * L - 1), dtype=np.float_) + np.random.randn(L, 2 * L - 1)

flm = ssht.forward(f, L, Reality=True, Method='GL')

f = ssht.inverse(flm, L, Reality=True, Method='GL')

flm = ssht.forward(f, L, Reality=True, Method='GL')

f_rec = ssht.inverse(flm, L, Reality=True, Method='GL')

if (np.mean(np.abs(f_rec - f)) < 1E-10):
    print('GL real transform test passed, max error:', np.max(np.abs(f_rec - f)))
else:
    print('GL real transform ***NOT*** test passed, max error:', np.max(
        np.abs(f_rec - f)))


def broken_function():
    raise Exception('This is broken')


class MyTestCase(unittest.TestCase):
    def test_forward_ndim(self):
        with self.assertRaises(ssht.ssht_input_error) as context:
            ssht.forward(flm, L)

        self.assertTrue('f must be 2D numpy array' in context.exception.args)

    def test_forward_method_type(self):
        with self.assertRaises(ssht.ssht_input_error) as context:
            ssht.forward(f, L, Method='DJ')

        self.assertTrue(
            'Method is not recognised, Methods are: MW, MW_pole, MWSS, DH and GL' in context.exception.args)

    def test_forward_spin_reality(self):
        with self.assertRaises(ssht.ssht_spin_error) as context:
            ssht.forward(f, L, Reality=True, Spin=2)

        self.assertTrue(
            'Reality set to True and Spin is not 0. However, spin signals must be complex.' in context.exception.args)

    def test_inverse_ndim(self):
        with self.assertRaises(ssht.ssht_input_error) as context:
            ssht.inverse(f, L)

        self.assertTrue('flm must be 1D numpy array' in context.exception.args)

    def test_inverse_method_type(self):
        with self.assertRaises(ssht.ssht_input_error) as context:
            ssht.inverse(flm, L, Method='DJ')

        self.assertTrue(
            'Method is not recognised, Methods are: MW, MW_pole, MWSS, DH and GL' in context.exception.args)

    def test_inverse_spin_reality(self):
        with self.assertRaises(ssht.ssht_spin_error) as context:
            ssht.forward(f, L, Reality=True, Spin=2)

        self.assertTrue(
            'Reality set to True and Spin is not 0. However, spin signals must be complex.' in context.exception.args)


if __name__ == "__main__":
    unittest.main()
