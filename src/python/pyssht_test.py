import numpy as np
from pyssht import *


L = 128
Spin = 0

flm = np.random.randn(L*L) + 1j*np.random.randn(L*L)

f = ssht_inverse(flm,L)

flm_rec = ssht_forward(f, L)

f_rec = ssht_inverse(flm_rec,L)

if (np.mean(np.abs(flm_rec-flm))<1E-14):
    print "MW complex transform test passed, max error: ", np.max(np.abs(flm_rec-flm))
else:
    print "MW complex transform ***NOT*** test passed, max error: ", np.max(np.abs(flm_rec-flm))
    
L = 128
Spin = 2

f = np.zeros((L,2*L-1), dtype=np.complex_) + np.random.randn(L,2*L-1)\
    + 1j*np.random.randn(L,2*L-1)

flm = ssht_forward(f, L, Spin=2)

f = ssht_inverse(flm,L,Spin=2)

flm = ssht_forward(f, L, Spin=2)

f_rec = ssht_inverse(flm,L,Spin=2)

if (np.mean(np.abs(f_rec-f))<1E-14):
    print "MW spin complex transform test passed, max error: ", np.max(np.abs(f_rec-f))
else:
    print "MW spin complex transform ***NOT*** test passed, max error: ", np.max(np.abs(f_rec-f))


L = 128
spin = 0

f = np.zeros((L,2*L-1), dtype=np.float_) + np.random.randn(L,2*L-1)

flm = ssht_forward(f, L,Reality=True)

f = ssht_inverse(flm,L,Reality=True)

flm = ssht_forward(f, L,Reality=True)

f_rec = ssht_inverse(flm,L,Reality=True)

if (np.mean(np.abs(f_rec-f))<1E-14):
    print "MW real transform test passed, max error: ", np.max(np.abs(f_rec-f))
else:
    print "MW real transform ***NOT*** test passed, max error: ", np.max(np.abs(f_rec-f))


L = 128
Spin = 0


flm = np.random.randn(L*L) + 1j*np.random.randn(L*L)

f, f_sp, phi_sp = ssht_inverse(flm,L, Method="MWSS")

flm_rec = ssht_forward((f, f_sp, phi_sp), L, Method="MWSS")

f_rec, f_sp_rec, phi_sp_rec  = ssht_inverse(flm,L, Method="MWSS")


if (np.mean(np.abs(flm_rec-flm))<1E-14):
    print "MWSS complex transform test passed, max error: ", np.max(np.abs(flm_rec-flm))
else:
    print "MWSS complex transform ***NOT*** test passed, max error: ", np.max(np.abs(flm_rec-flm))


    
L = 128
Spin = 2


flm = np.random.randn(L*L) + 1j*np.random.randn(L*L)
flm[0:4] = 0.0 + 0.0j

f, f_sp, phi_sp = ssht_inverse_mw_complex_pole(flm,L, Spin)

flm_rec = ssht_forward_mw_complex_pole(f, f_sp, phi_sp, L, Spin)

f_rec, f_sp_rec, phi_sp_rec  = ssht_inverse_mw_complex_pole(flm,L, Spin)

if (np.mean(np.abs(flm_rec-flm))<1E-14 and np.abs(f_sp_rec-f_sp)<1E-12 and np.abs(phi_sp_rec-phi_sp)<1E-12):
    print "MWSS spin complex transform test passed, max error: ", np.max(np.abs(flm_rec-flm))
else:
    print "MWSS spin complex transform ***NOT*** test passed, max error: ", np.max(np.abs(flm_rec-flm))

    
    
L = 128
Spin = 0

f = np.random.randn(L-1,2*L)

f_sp = np.random.randn()

flm = ssht_forward_mw_real_pole(f, f_sp, L)

f, f_sp = ssht_inverse_mw_real_pole(flm,L)

flm = ssht_forward_mw_real_pole(f, f_sp, L)

f_rec, f_sp_rec  = ssht_inverse_mw_real_pole(flm,L)

f_dif = f-f_rec

#print f_sp, f_sp_rec, np.abs(f_sp_rec-f_sp)
#print f_rec[L-1,0:10]
#print f[L-1,0:10]
#print f_rec[0,0:10]
#print f[0,0:10]
#print np.abs(f_dif[L-1,:])

if (np.mean(np.abs(f_rec-f))<1E-14 and np.abs(f_sp_rec-f_sp)<1E-12):
    print "MWSS real transform test passed, max error: ", np.max(np.abs(f_rec-f))
else:
    print "MWSS real transform ***NOT*** test passed, max error: ", np.max(np.abs(f_rec-f))

    
    
L = 128
Spin = 0

f = np.zeros((2*L,2*L-1), dtype=np.complex_) + np.random.randn(2*L,2*L-1)\
    + 1j*np.random.randn(2*L,2*L-1)

flm = ssht_forward(f, L, Method="DH")

f = ssht_inverse(flm,L, Method="DH")

flm = ssht_forward(f, L, Method="DH")

f_rec = ssht_inverse(flm,L, Method="DH")

if (np.mean(np.abs(f_rec-f))<1E-10):
    print "DH complex transform test passed, max error: ", np.max(np.abs(f_rec-f))
else:
    print "DH complex transform ***NOT*** test passed, max error: ", np.max(np.abs(f_rec-f))
    

L = 128
Spin = 2

f = np.zeros((2*L,2*L-1), dtype=np.complex_) + np.random.randn(2*L,2*L-1)\
    + 1j*np.random.randn(2*L,2*L-1)

flm = ssht_forward(f, L, Spin=Spin, Method="DH")

f = ssht_inverse(flm,L, Spin=Spin, Method="DH")

flm = ssht_forward(f, L, Spin=Spin, Method="DH")

f_rec = ssht_inverse(flm,L, Spin=Spin, Method="DH")

if (np.mean(np.abs(f_rec-f))<1E-10):
    print "DH spin complex transform test passed, max error: ", np.max(np.abs(f_rec-f))
else:
    print "DH spin complex transform ***NOT*** test passed, max error: ", np.max(np.abs(f_rec-f))


L = 128
spin = 0

f = np.zeros((2*L,2*L-1), dtype=np.float_) + np.random.randn(2*L,2*L-1)

flm = ssht_forward(f, L, Method="DH", Reality=True)

f = ssht_inverse(flm,L, Method="DH", Reality=True)

flm = ssht_forward(f, L, Method="DH", Reality=True)

f_rec = ssht_inverse(flm,L, Method="DH", Reality=True)

if (np.mean(np.abs(f_rec-f))<1E-10):
    print "DH real transform test passed, max error: ", np.max(np.abs(f_rec-f))
else:
    print "DH real transform ***NOT*** test passed, max error: ", np.max(np.abs(f_rec-f))


 
L = 128
Spin = 0

f = np.zeros((L,2*L-1), dtype=np.complex_) + np.random.randn(L,2*L-1)\
    + 1j*np.random.randn(L,2*L-1)

flm = ssht_forward(f, L, Method="GL")

f = ssht_inverse(flm,L,Method="GL")

flm = ssht_forward(f, L, Method="GL")

f_rec = ssht_inverse(flm,L,Method="GL")

if (np.mean(np.abs(f_rec-f))<1E-10):
    print "GL complex transform test passed, max error: ", np.max(np.abs(f_rec-f))
else:
    print "GL complex transform ***NOT*** test passed, max error: ", np.max(np.abs(f_rec-f))


L = 128
Spin = 2

f = np.zeros((L,2*L-1), dtype=np.complex_) + np.random.randn(L,2*L-1)\
    + 1j*np.random.randn(L,2*L-1)

flm = ssht_forward(f, L, Spin=Spin, Method="GL")

f = ssht_inverse(flm,L,Spin=Spin, Method="GL")

flm = ssht_forward(f, L, Spin=Spin, Method="GL")

f_rec = ssht_inverse(flm,L,Spin=Spin, Method="GL")

if (np.mean(np.abs(f_rec-f))<1E-10):
    print "GL spin complex transform test passed, max error: ", np.max(np.abs(f_rec-f))
else:
    print "GL spin complex transform ***NOT*** test passed, max error: ", np.max(np.abs(f_rec-f))


L = 128
spin = 0

f = np.zeros((L,2*L-1), dtype=np.float_) + np.random.randn(L,2*L-1)

flm = ssht_forward(f, L, Reality=True, Method="GL")

f = ssht_inverse(flm,L, Reality=True, Method="GL")

flm = ssht_forward(f, L, Reality=True, Method="GL")

f_rec = ssht_inverse(flm,L, Reality=True, Method="GL")

if (np.mean(np.abs(f_rec-f))<1E-10):
    print "GL real transform test passed, max error: ", np.max(np.abs(f_rec-f))
else:
    print "GL real transform ***NOT*** test passed, max error: ", np.max(np.abs(f_rec-f))

     
import unittest

def broken_function():
    raise Exception('This is broken')

class MyTestCase(unittest.TestCase):
    
    def test_forward_ndim(self):
        with self.assertRaises(ssht_input_error) as context:
            ssht_forward(flm,L)

        self.assertTrue('f must be 2D numpy array' in context.exception)

    def test_forward_method_type(self):
        with self.assertRaises(ssht_input_error) as context:
            ssht_forward(f,L,Method="DJ")

        self.assertTrue('Method is not recognised, Methods are: MW, MWSS, DH and GL' in context.exception)

    def test_forward_spin_reality(self):
        with self.assertRaises(ssht_spin_error) as context:
            ssht_forward(f,L,Reality=True,Spin=2)

        self.assertTrue('Reality set to True and Spin is not 0. However spin signals must be complex.'
                        in context.exception)
        
        
    def test_inverse_ndim(self):
        with self.assertRaises(ssht_input_error) as context:
            ssht_inverse(f,L)

        self.assertTrue('flm must be 1D numpy array' in context.exception)

    def test_inverse_method_type(self):
        with self.assertRaises(ssht_input_error) as context:
            ssht_inverse(flm,L,Method="DJ")

        self.assertTrue('Method is not recognised, Methods are: MW, MWSS, DH and GL' in context.exception)
        
    def test_inverse_spin_reality(self):
        with self.assertRaises(ssht_spin_error) as context:
            ssht_forward(f,L,Reality=True,Spin=2)

        self.assertTrue('Reality set to True and Spin is not 0. However spin signals must be complex.'
                        in context.exception)
        

unittest.main()


