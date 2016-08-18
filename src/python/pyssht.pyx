# import both numpy and the Cython declarations for numpy
import numpy as np
cimport numpy as np

#----------------------------------------------------------------------------------------------------#

cdef extern from "ssht.h":

        double ssht_sampling_mw_t2theta(int t, int L); # add the other samplings, ploting, c2s etc; dl, rotations, adjoints
        double ssht_sampling_mw_p2phi(int p, int L);
        int ssht_sampling_mw_n(int L);
        int ssht_sampling_mw_ntheta(int L);
        int ssht_sampling_mw_nphi(int L);
        double* ssht_dl_calloc(int L, ssht_dl_size_t dl_size);
        ctypedef enum ssht_dl_size_t:
                SSHT_DL_QUARTER_EXTENDED, SSHT_DL_HALF, SSHT_DL_FULL

        void ssht_dl_beta_risbo_full_table(double *dl, double beta, int L,
                                           ssht_dl_size_t dl_size,
                                           int el, double *sqrt_tbl);

        void ssht_dl_beta_risbo_half_table(double *dl, double beta, int L,
                                           ssht_dl_size_t dl_size,
                                           int el, double *sqrt_tbl, double *signs);

# I included
        ctypedef enum  ssht_dl_method_t:
                SSHT_DL_RISBO, SSHT_DL_TRAPANI
        void ssht_core_mw_forward_sov_conv_sym(double complex *flm, const double complex *f,
				       int L, int spin,
				       ssht_dl_method_t dl_method,
				       int verbosity);
        void ssht_core_mw_inverse_sov_sym(double complex *f, const double complex *flm,
                                  int L, int spin,
                                  ssht_dl_method_t dl_method,
                                  int verbosity);
        void ssht_core_mw_inverse_sov_sym_real(double *f, const double complex *flm,
                                  int L,
                                  ssht_dl_method_t dl_method,
                                  int verbosity);
        void ssht_core_mw_forward_sov_conv_sym_real(double complex *flm, const double *f,
                                  int L,
                                  ssht_dl_method_t dl_method,
                                  int verbosity);
        void ssht_core_mw_inverse_sov_sym_pole(double complex *f,
                                  double complex *f_sp, double *phi_sp,
                                  const double complex *flm,
                                  int L, int spin,
                                  ssht_dl_method_t dl_method,
                                  int verbosity);
        void ssht_core_mw_inverse_sov_sym_real_pole(double *f,
                                  double *f_sp,
                                  const double complex *flm,
                                  int L,
                                  ssht_dl_method_t dl_method,
                                  int verbosity);
        void ssht_core_mw_forward_sov_conv_sym_pole(double complex *flm, const double complex *f, # change to ss 
                                  double complex f_sp, double phi_sp,
                                  int L, int spin,
                                  ssht_dl_method_t dl_method,
                                  int verbosity);
        void ssht_core_mw_forward_sov_conv_sym_real_pole(double complex *flm,
                                  const double *f,
                                  double f_sp,
                                  int L,
                                  ssht_dl_method_t dl_method,
                                  int verbosity);
        void ssht_core_mw_inverse_sov_sym_ss(double complex *f, const double complex *flm,
              int L, int spin,
              ssht_dl_method_t dl_method,
              int verbosity);
        void ssht_core_mw_inverse_sov_sym_ss_real(double *f, const double complex *flm,
              int L,
              ssht_dl_method_t dl_method,
              int verbosity);
        void ssht_core_mw_forward_sov_conv_sym_ss(double complex *flm, const double complex *f,
              int L, int spin,
              ssht_dl_method_t dl_method,
              int verbosity);
        void ssht_core_mw_forward_sov_conv_sym_ss_real(double complex *flm, const double *f,
              int L,
              ssht_dl_method_t dl_method,
              int verbosity);
        void ssht_core_dh_inverse_sov(double complex *f, const double complex *flm,
                                  int L, int spin, int verbosity);
        void ssht_core_dh_inverse_sov_real(double *f, const double complex *flm,
                                  int L, int verbosity);
        void ssht_core_dh_forward_sov(double complex *flm, const double complex *f,
                                  int L, int spin, int verbosity);
        void ssht_core_dh_forward_sov_real(double complex *flm, const double *f,
                                  int L, int verbosity);
        void ssht_core_gl_inverse_sov(double complex *f, const double complex *flm,
                                  int L, int spin, int verbosity);
        void ssht_core_gl_inverse_sov_real(double *f, const double complex *flm,
                                  int L, int verbosity);
        void ssht_core_gl_forward_sov(double complex *flm, const double complex *f,
                                  int L, int spin, int verbosity);
        void ssht_core_gl_forward_sov_real(double complex *flm, const double *f,
                                  int L, int verbosity);
        double ssht_sampling_mw_t2theta(int t, int L);
        double ssht_sampling_mw_p2phi(int p, int L);
        double ssht_sampling_mw_ss_t2theta(int t, int L);
        double ssht_sampling_mw_ss_p2phi(int p, int L);
        double ssht_sampling_dh_t2theta(int t, int L);
        double ssht_sampling_dh_p2phi(int p, int L);
        void ssht_sampling_gl_thetas_weights(double *thetas, double *weights, int L);
        double ssht_sampling_gl_p2phi(int p, int L);



#----------------------------------------------------------------------------------------------------#

def ssht_forward_mw_complex(np.ndarray[ double complex, ndim=2, mode="c"] f_mw_c not None,int L,int spin):

        cdef ssht_dl_method_t dl_method = SSHT_DL_RISBO;
        f_lm = np.empty([L * L,], dtype=complex)
        ssht_core_mw_forward_sov_conv_sym(<double complex*> np.PyArray_DATA(f_lm),<const double complex*> np.PyArray_DATA(f_mw_c), L, spin, dl_method, 0);
        return f_lm
        
def ssht_inverse_mw_complex(np.ndarray[ double complex, ndim=1, mode="c"] f_lm not None, int L, int spin):

        cdef ssht_dl_method_t dl_method = SSHT_DL_RISBO;
        f_mw_c = np.empty([L,2*L-1,], dtype=complex)
        ssht_core_mw_inverse_sov_sym(<double complex*> np.PyArray_DATA(f_mw_c),<const double complex*> np.PyArray_DATA(f_lm), L, spin, dl_method, 0);
        return f_mw_c

#----------------------------------------------------------------------------------------------------#

def ssht_forward_mw_real(np.ndarray[ double, ndim=2, mode="c"] f_mw_r not None,int L):

        cdef ssht_dl_method_t dl_method = SSHT_DL_RISBO;
        f_lm = np.empty([L * L,], dtype=complex)
        ssht_core_mw_forward_sov_conv_sym_real(<double complex*> np.PyArray_DATA(f_lm),<const double*> np.PyArray_DATA(f_mw_r), L, dl_method, 0);
        return f_lm
        
def ssht_inverse_mw_real(np.ndarray[ double complex, ndim=1, mode="c"] f_lm not None, int L):

        cdef ssht_dl_method_t dl_method = SSHT_DL_RISBO;
        f_mw_r = np.empty([L,2*L-1,], dtype=np.float_)
        ssht_core_mw_inverse_sov_sym_real(<double*> np.PyArray_DATA(f_mw_r),<const double complex*> np.PyArray_DATA(f_lm), L, dl_method, 0);
        return f_mw_r

#----------------------------------------------------------------------------------------------------#

def ssht_forward_mwss_complex(np.ndarray[ double complex, ndim=2, mode="c"] f_mwss_c not None,int L,int spin):

        cdef ssht_dl_method_t dl_method = SSHT_DL_RISBO;
        f_lm = np.empty([L * L,], dtype=complex)
        ssht_core_mw_forward_sov_conv_sym_ss(<double complex*> np.PyArray_DATA(f_lm),<const double complex*> np.PyArray_DATA(f_mwss_c), L, spin, dl_method, 0);
        return f_lm
        
def ssht_inverse_mwss_complex(np.ndarray[ double complex, ndim=1, mode="c"] f_lm not None, int L, int spin):

        cdef ssht_dl_method_t dl_method = SSHT_DL_RISBO;
        f_mwss_c = np.empty([L+1,2*L,], dtype=complex)
        ssht_core_mw_inverse_sov_sym_ss(<double complex*> np.PyArray_DATA(f_mwss_c),<const double complex*> np.PyArray_DATA(f_lm), L, spin, dl_method, 0);
        return f_mwss_c

#----------------------------------------------------------------------------------------------------#

def ssht_forward_mwss_real(np.ndarray[ double, ndim=2, mode="c"] f_mwss_r not None,int L):

        cdef ssht_dl_method_t dl_method = SSHT_DL_RISBO;
        f_lm = np.empty([L * L,], dtype=complex)
        ssht_core_mw_forward_sov_conv_sym_ss_real(<double complex*> np.PyArray_DATA(f_lm),<const double*> np.PyArray_DATA(f_mwss_r), L, dl_method, 0);
        return f_lm
        
def ssht_inverse_mwss_real(np.ndarray[ double complex, ndim=1, mode="c"] f_lm not None, int L):

        cdef ssht_dl_method_t dl_method = SSHT_DL_RISBO;
        f_mwss_r = np.empty([L+1,2*L,], dtype=np.float_)
        ssht_core_mw_inverse_sov_sym_ss_real(<double*> np.PyArray_DATA(f_mwss_r),<const double complex*> np.PyArray_DATA(f_lm), L, dl_method, 0);
        return f_mwss_r

#----------------------------------------------------------------------------------------------------#

def ssht_forward_mw_complex_pole(np.ndarray[ double complex, ndim=2, mode="c"] f_mw_c not None, double complex f_sp, double phi_sp, int L,int spin):

        cdef ssht_dl_method_t dl_method = SSHT_DL_RISBO;
        f_lm = np.empty([L * L,], dtype=complex)
        ssht_core_mw_forward_sov_conv_sym_pole(<double complex*> np.PyArray_DATA(f_lm),<const double complex*> np.PyArray_DATA(f_mw_c),  f_sp,  phi_sp, L, spin, dl_method, 0);
        return f_lm
        
def ssht_inverse_mw_complex_pole(np.ndarray[ double complex, ndim=1, mode="c"] f_lm not None, int L, int spin):

        cdef ssht_dl_method_t dl_method = SSHT_DL_RISBO;
        f_mw_c = np.empty([L-1,(2*L-1),], dtype=complex)
        cdef double complex f_sp
        cdef double phi_sp
        ssht_core_mw_inverse_sov_sym_pole(<double complex*> np.PyArray_DATA(f_mw_c),  &f_sp,  &phi_sp, <const double complex*> np.PyArray_DATA(f_lm), L, spin, dl_method, 0);
        return f_mw_c, f_sp, phi_sp

#----------------------------------------------------------------------------------------------------#

def ssht_forward_mw_real_pole(np.ndarray[ double, ndim=2, mode="c"] f_mw_r not None, double f_sp,  int L):

        cdef ssht_dl_method_t dl_method = SSHT_DL_RISBO;
        f_lm = np.empty([L * L,], dtype=complex)
        ssht_core_mw_forward_sov_conv_sym_real_pole(<double complex*> np.PyArray_DATA(f_lm),<const double*> np.PyArray_DATA(f_mw_r), f_sp, L, dl_method, 0);
        return f_lm
        
def ssht_inverse_mw_real_pole(np.ndarray[ double complex, ndim=1, mode="c"] f_lm not None, int L):

        cdef ssht_dl_method_t dl_method = SSHT_DL_RISBO;
        f_mw_r = np.empty([L-1,2*L-1,], dtype=np.float_)
        cdef double f_sp
        ssht_core_mw_inverse_sov_sym_real_pole(<double*> np.PyArray_DATA(f_mw_r), &f_sp, <const double complex*> np.PyArray_DATA(f_lm), L, dl_method, 0);
        return f_mw_r, f_sp

#----------------------------------------------------------------------------------------------------#

def ssht_forward_dh_complex(np.ndarray[ double complex, ndim=2, mode="c"] f_dh_c not None,int L,int spin):

        f_lm = np.empty([L * L,], dtype=complex)
        ssht_core_dh_forward_sov(<double complex*> np.PyArray_DATA(f_lm),<const double complex*> np.PyArray_DATA(f_dh_c), L, spin, 0);
        return f_lm
        
def ssht_inverse_dh_complex(np.ndarray[ double complex, ndim=1, mode="c"] f_lm not None, int L, int spin):

        f_dh_c = np.empty([2*L,2*L-1,], dtype=complex)
        ssht_core_dh_inverse_sov(<double complex*> np.PyArray_DATA(f_dh_c),<const double complex*> np.PyArray_DATA(f_lm), L, spin, 0);
        return f_dh_c

#----------------------------------------------------------------------------------------------------#

def ssht_forward_dh_real(np.ndarray[ double, ndim=2, mode="c"] f_dh_r not None,int L):

        f_lm = np.empty([L * L,], dtype=complex)
        ssht_core_dh_forward_sov_real(<double complex*> np.PyArray_DATA(f_lm),<const double*> np.PyArray_DATA(f_dh_r), L,  0);
        return f_lm
        
def ssht_inverse_dh_real(np.ndarray[ double complex, ndim=1, mode="c"] f_lm not None, int L):

        f_dh_r = np.empty([2*L,2*L-1,], dtype=np.float_)
        ssht_core_dh_inverse_sov_real(<double*> np.PyArray_DATA(f_dh_r),<const double complex*> np.PyArray_DATA(f_lm), L, 0);
        return f_dh_r


#----------------------------------------------------------------------------------------------------#

def ssht_forward_gl_complex(np.ndarray[ double complex, ndim=2, mode="c"] f_gl_c not None,int L,int spin):

        cdef ssht_dl_method_t dl_method = SSHT_DL_RISBO;
        f_lm = np.empty([L * L,], dtype=complex)
        ssht_core_gl_forward_sov(<double complex*> np.PyArray_DATA(f_lm),<const double complex*> np.PyArray_DATA(f_gl_c), L, spin, 0);
        return f_lm
        
def ssht_inverse_gl_complex(np.ndarray[ double complex, ndim=1, mode="c"] f_lm not None, int L, int spin):

        cdef ssht_dl_method_t dl_method = SSHT_DL_RISBO;
        f_gl_c = np.empty([L,2*L-1,], dtype=complex)
        ssht_core_gl_inverse_sov(<double complex*> np.PyArray_DATA(f_gl_c),<const double complex*> np.PyArray_DATA(f_lm), L, spin, 0);
        return f_gl_c

#----------------------------------------------------------------------------------------------------#

def ssht_forward_gl_real(np.ndarray[ double, ndim=2, mode="c"] f_gl_r not None,int L):

        cdef ssht_dl_method_t dl_method = SSHT_DL_RISBO;
        f_lm = np.empty([L * L,], dtype=complex)
        ssht_core_gl_forward_sov_real(<double complex*> np.PyArray_DATA(f_lm),<const double*> np.PyArray_DATA(f_gl_r), L, 0);
        return f_lm
        
def ssht_inverse_gl_real(np.ndarray[ double complex, ndim=1, mode="c"] f_lm not None, int L):

        cdef ssht_dl_method_t dl_method = SSHT_DL_RISBO;
        f_gl_r = np.empty([L,2*L-1,], dtype=np.float_)
        ssht_core_gl_inverse_sov_real(<double*> np.PyArray_DATA(f_gl_r),<const double complex*> np.PyArray_DATA(f_lm), L, 0);
        return f_gl_r



#----------------------------------------------------------------------------------------------------#

# some error handling

class ssht_input_error(TypeError):
    '''Raise this when inputs are the wrong type'''

class ssht_spin_error(ValueError):
    '''Raise this when the spin is none zero but Reality is true'''


#----------------------------------------------------------------------------------------------------#

# easy to use functions to perform forward and backward transfroms

def ssht_forward(f, int L, Spin=0, Method='MW', Reality=False):
    # Checks

    if Method == 'MW_pole':
        if Reality:
            f, f_sp = f
        else:
            f, f_sp, phi_sp = f

    if f.ndim != 2:
      raise ssht_input_error('f must be 2D numpy array')

    if not(Method == 'MW' or Method == 'MW_pole' or Method == 'MWSS' or Method == 'DH' or Method == 'GL'):
        raise ssht_input_error('Method is not recognised, Methods are: MW, MW_pole, MWSS, DH and GL')

    if Spin != 0 and Reality == True :
        raise ssht_spin_error('Reality set to True and Spin is not 0. However, spin signals must be complex.')

    if f.dtype == np.float_ and Reality == False:
        print 'Real signal given but Reality flag is False. Set Reality = True to improve performance'
        f_new = np.empty((L,2*L-1), dtype=np.complex_)
        f_new = f + 1j*np.zeros((L,2*L-1), dtype=np.float_)
        f = f_new

    if f.dtype == np.complex_ and Reality == True:
        print 'Complex signal given but Reality flag is True. Ignoring complex component'
        f_new = np.real(f)
        f = f_new.copy(order='c')
        
    # do correct transform
    if Method == 'MW':
        if Reality:
            flm = ssht_forward_mw_real(f,L)
        else:
            flm = ssht_forward_mw_complex(f,L,Spin)
            
    if Method == 'MW_pole':
        if Reality:
            flm = ssht_forward_mw_real_pole(f,f_sp,L)
        else:
            flm = ssht_forward_mw_complex_pole(f,f_sp,phi_sp,L,Spin)
            
    if Method == 'MWSS':
        if Reality:
            flm = ssht_forward_mwss_real(f,L)
        else:
            flm = ssht_forward_mwss_complex(f,L,Spin)

    if Method == 'DH':
        if Reality:
            flm = ssht_forward_dh_real(f,L)
        else:
            flm = ssht_forward_dh_complex(f,L,Spin)
            
    if Method == 'GL':
        if Reality:
            flm = ssht_forward_gl_real(f,L)
        else:
            flm = ssht_forward_gl_complex(f,L,Spin)
            
            
    return flm


def ssht_inverse(flm, L, Spin=0, Method='MW', Reality=False):
    # Checks
    if flm.ndim != 1:
        raise ssht_input_error('flm must be 1D numpy array')

    if not(Method == 'MW' or Method == 'MW_pole' or Method == 'MWSS' or Method == 'DH' or Method == "GL"):
        raise ssht_input_error('Method is not recognised, Methods are: MW, MW_pole, MWSS, DH and GL')
    
    if Spin != 0 and Reality == True :
        raise ssht_spin_error('Reality set to True and Spin is not 0. However, spin signals must be complex.')
    
    # do correct transform
    if Method == 'MW':
        if Reality:
            f = ssht_inverse_mw_real(flm,L)
        else:
            f = ssht_inverse_mw_complex(flm,L,Spin)
            

    if Method == 'MW_pole':
        if Reality:
            f = ssht_inverse_mw_real_pole(flm,L)
        else:
            f = ssht_inverse_mw_complex_pole(flm,L,Spin)

    if Method == 'MWSS':
        if Reality:
            f = ssht_inverse_mwss_real(flm,L)
        else:
            f = ssht_inverse_mwss_complex(flm,L,Spin)
            
    if Method == 'DH':
        if Reality:
            f = ssht_inverse_dh_real(flm,L)
        else:
            f = ssht_inverse_dh_complex(flm,L,Spin)
            
    if Method == 'GL':
        if Reality:
            f = ssht_inverse_gl_real(flm,L)
        else:
            f = ssht_inverse_gl_complex(flm,L,Spin)
            
            
            
    return f

#----------------------------------------------------------------------------------------------------#

# index to ell em and back function

def isqrt(int n):
    cdef int square = 1, delta = 3
    while square < n:
        square += delta
        delta  += 2
    return (delta/2 -1)

def ssht_elm2ind( int el, int m):

  return el * el + el + m


def ssht_ind2elm(int ind):

  cdef int ell, em
  el = isqrt(ind)
  em = ind - (el)*(el) - (el);

  return el, em

# get the shape of the signal on the sphere for different sampling theorems

def ssht_sample_shape(int L, Method = 'MW'):
  if not(Method == 'MW' or Method == 'MW_pole' or Method == 'MWSS' or Method == 'DH' or Method == "GL"):
      raise ssht_input_error('Method is not recognised, Methods are: MW, MW_pole, MWSS, DH and GL')

  cdef int n_theta, n_phi
  if Method == 'MW':
    n_theta = L
    n_phi   = 2*L-1

  if Method == 'MW_pole':
    n_theta = L-1
    n_phi   = 2*L-1

  if Method == 'MWSS':
    n_theta = L+1
    n_phi   = 2*L

  if Method == 'GL':
    n_theta = L
    n_phi   = 2*L-1

  if Method == 'DH':
    n_theta = 2*L
    n_phi   = 2*L-1



  return (n_theta, n_phi)

def ssht_sample_positions(int L, Method = 'MW', Grid=False):
  if not(Method == 'MW' or Method == 'MW_pole' or Method == 'MWSS' or Method == 'DH' or Method == "GL"):
    raise ssht_input_error('Method is not recognised, Methods are: MW, MW_pole, MWSS, DH and GL')

  n_theta, n_phi = ssht_sample_shape(L, Method=Method)
  thetas = np.empty(n_theta, dtype=np.float_)
  phis   = np.empty(n_phi,   dtype=np.float_)


  if Method == 'MW':
    for i in range(n_theta):
      thetas[i] = ssht_sampling_mw_t2theta(i, L)
    for i in range(n_phi):
      phis[i] = ssht_sampling_mw_p2phi(i, L)

  if Method == 'MWSS':
    for i in range(n_theta):
      thetas[i] = ssht_sampling_mw_ss_t2theta(i, L)
    for i in range(n_phi):
      phis[i] = ssht_sampling_mw_ss_p2phi(i, L)

  if Method == 'DH':
    for i in range(n_theta):
      thetas[i] = ssht_sampling_dh_t2theta(i, L)
    for i in range(n_phi):
      phis[i] = ssht_sampling_dh_p2phi(i, L)

  if Method == 'GL':
    weights_unused = np.empty(L, dtype=np.float_)
    ssht_sampling_gl_thetas_weights(<double*> np.PyArray_DATA(thetas), <double*> np.PyArray_DATA(weights_unused), L)
    for i in range(n_phi):
      phis[i] = ssht_sampling_gl_p2phi(i, L)


  if Grid:
    phis, thetas = np.meshgrid(phis, thetas)

  return thetas, phis






