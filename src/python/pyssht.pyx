# import both numpy and the Cython declarations for numpy
import numpy as np
cimport numpy as np

import matplotlib.pyplot as plt
from matplotlib import cm, colors, colorbar, gridspec
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.axes_grid1 import make_axes_locatable

#----------------------------------------------------------------------------------------------------#

cdef enum METHOD_TYPE:
  MW, MWSS, DH, GL

#----------------------------------------------------------------------------------------------------#

cdef extern from "ssht.h":

        double ssht_sampling_mw_t2theta(int t, int L); #  ; adjoints
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

# adjoints

        void ssht_adjoint_mw_inverse_sov_sym(double complex *flm, 
             double complex *f, 
             int L, int spin, 
             ssht_dl_method_t dl_method,
             int verbosity);
        void ssht_adjoint_mw_inverse_sov_sym_real(double complex *flm, 
            double *f, 
            int L,
            ssht_dl_method_t dl_method, 
            int verbosity);
        void ssht_adjoint_mw_forward_sov_sym(double complex *f, 
             double complex *flm,
             int L, int spin,
             ssht_dl_method_t dl_method,
             int verbosity);
        void ssht_adjoint_mw_forward_sov_sym_real(double *f, 
            double complex *flm,
            int L,
            ssht_dl_method_t dl_method,
            int verbosity);

        void ssht_adjoint_mw_inverse_sov_sym_pole(double complex *flm, double complex *f,
            double complex f_sp, double phi_sp,
            int L, int spin, 
            ssht_dl_method_t dl_method,
            int verbosity);
        void ssht_adjoint_mw_inverse_sov_sym_real_pole(double complex *flm, 
                 double *f, 
                 double f_sp,
                 int L, 
                 ssht_dl_method_t dl_method,
                 int verbosity);
        void ssht_adjoint_mw_forward_sov_sym_pole(double complex *f, 
            double complex *f_sp, double *phi_sp,
            double complex *flm, 
            int L, int spin, 
            ssht_dl_method_t dl_method,
            int verbosity);
        void ssht_adjoint_mw_forward_sov_sym_real_pole(double *f, 
                 double *f_sp,
                 double complex *flm, 
                 int L, 
                 ssht_dl_method_t dl_method, 
                 int verbosity);
        void ssht_adjoint_mw_inverse_sov_sym_ss(double complex *flm, double complex *f, 
          int L, int spin, 
          ssht_dl_method_t dl_method,
          int verbosity);
        void ssht_adjoint_mw_inverse_sov_sym_ss_real(double complex *flm, double *f, 
               int L, 
               ssht_dl_method_t dl_method, 
               int verbosity);
        void ssht_adjoint_mw_forward_sov_sym_ss(double complex *f, double complex *flm,
          int L, int spin,
          ssht_dl_method_t dl_method,
          int verbosity);
        void ssht_adjoint_mw_forward_sov_sym_ss_real(double *f, 
               double complex *flm,
               int L,
               ssht_dl_method_t dl_method,
               int verbosity);



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


def ssht_inverse_mw_complex_adjoint(np.ndarray[ double complex, ndim=2, mode="c"] f_mw_c not None,int L,int spin):

        cdef ssht_dl_method_t dl_method = SSHT_DL_RISBO;
        f_lm = np.empty([L * L,], dtype=complex)
        ssht_adjoint_mw_inverse_sov_sym(<double complex*> np.PyArray_DATA(f_lm),<const double complex*> np.PyArray_DATA(f_mw_c), L, spin, dl_method, 0);
        return f_lm
        
def ssht_forward_mw_complex_adjoint(np.ndarray[ double complex, ndim=1, mode="c"] f_lm not None, int L, int spin):

        cdef ssht_dl_method_t dl_method = SSHT_DL_RISBO;
        f_mw_c = np.empty([L,2*L-1,], dtype=complex)
        ssht_adjoint_mw_forward_sov_sym(<double complex*> np.PyArray_DATA(f_mw_c),<const double complex*> np.PyArray_DATA(f_lm), L, spin, dl_method, 0);
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

def ssht_inverse_mw_real_adjoint(np.ndarray[ double, ndim=2, mode="c"] f_mw_r not None,int L):

        cdef ssht_dl_method_t dl_method = SSHT_DL_RISBO;
        f_lm = np.empty([L * L,], dtype=complex)
        ssht_adjoint_mw_inverse_sov_sym_real(<double complex*> np.PyArray_DATA(f_lm),<const double*> np.PyArray_DATA(f_mw_r), L, dl_method, 0);
        return f_lm
        
def ssht_forward_mw_real_adjoint(np.ndarray[ double complex, ndim=1, mode="c"] f_lm not None, int L):

        cdef ssht_dl_method_t dl_method = SSHT_DL_RISBO;
        f_mw_r = np.empty([L,2*L-1,], dtype=np.float_)
        ssht_adjoint_mw_forward_sov_sym_real(<double*> np.PyArray_DATA(f_mw_r),<const double complex*> np.PyArray_DATA(f_lm), L, dl_method, 0);
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


def ssht_inverse_mwss_complex_adjoint(np.ndarray[ double complex, ndim=2, mode="c"] f_mwss_c not None,int L,int spin):

        cdef ssht_dl_method_t dl_method = SSHT_DL_RISBO;
        f_lm = np.empty([L * L,], dtype=complex)
        ssht_adjoint_mw_inverse_sov_sym_ss(<double complex*> np.PyArray_DATA(f_lm),<const double complex*> np.PyArray_DATA(f_mwss_c), L, spin, dl_method, 0);
        return f_lm
        
def ssht_forward_mwss_complex_adjoint(np.ndarray[ double complex, ndim=1, mode="c"] f_lm not None, int L, int spin):

        cdef ssht_dl_method_t dl_method = SSHT_DL_RISBO;
        f_mwss_c = np.empty([L+1,2*L,], dtype=complex)
        ssht_adjoint_mw_forward_sov_sym_ss(<double complex*> np.PyArray_DATA(f_mwss_c),<const double complex*> np.PyArray_DATA(f_lm), L, spin, dl_method, 0);
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

def ssht_inverse_mwss_real_adjoint(np.ndarray[ double, ndim=2, mode="c"] f_mwss_r not None,int L):

        cdef ssht_dl_method_t dl_method = SSHT_DL_RISBO;
        f_lm = np.empty([L * L,], dtype=complex)
        ssht_adjoint_mw_inverse_sov_sym_ss_real(<double complex*> np.PyArray_DATA(f_lm),<const double*> np.PyArray_DATA(f_mwss_r), L, dl_method, 0);
        return f_lm
        
def ssht_forward_mwss_real_adjoint(np.ndarray[ double complex, ndim=1, mode="c"] f_lm not None, int L):

        cdef ssht_dl_method_t dl_method = SSHT_DL_RISBO;
        f_mwss_r = np.empty([L+1,2*L,], dtype=np.float_)
        ssht_adjoint_mw_forward_sov_sym_ss_real(<double*> np.PyArray_DATA(f_mwss_r),<const double complex*> np.PyArray_DATA(f_lm), L, dl_method, 0);
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


def ssht_inverse_mw_complex_pole_adjoint(np.ndarray[ double complex, ndim=2, mode="c"] f_mw_c not None, double complex f_sp, double phi_sp, int L,int spin):

        cdef ssht_dl_method_t dl_method = SSHT_DL_RISBO;
        f_lm = np.empty([L * L,], dtype=complex)
        ssht_adjoint_mw_inverse_sov_sym_pole(<double complex*> np.PyArray_DATA(f_lm),<const double complex*> np.PyArray_DATA(f_mw_c),  f_sp,  phi_sp, L, spin, dl_method, 0);
        return f_lm
        
def ssht_forward_mw_complex_pole_adjoint(np.ndarray[ double complex, ndim=1, mode="c"] f_lm not None, int L, int spin):

        cdef ssht_dl_method_t dl_method = SSHT_DL_RISBO;
        f_mw_c = np.empty([L-1,(2*L-1),], dtype=complex)
        cdef double complex f_sp
        cdef double phi_sp
        ssht_adjoint_mw_forward_sov_sym_pole(<double complex*> np.PyArray_DATA(f_mw_c),  &f_sp,  &phi_sp, <const double complex*> np.PyArray_DATA(f_lm), L, spin, dl_method, 0);
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

def ssht_inverse_mw_real_pole_adjoint(np.ndarray[ double, ndim=2, mode="c"] f_mw_r not None, double f_sp,  int L):

        cdef ssht_dl_method_t dl_method = SSHT_DL_RISBO;
        f_lm = np.empty([L * L,], dtype=complex)
        ssht_adjoint_mw_inverse_sov_sym_real_pole(<double complex*> np.PyArray_DATA(f_lm),<const double*> np.PyArray_DATA(f_mw_r), f_sp, L, dl_method, 0);
        return f_lm
        
def ssht_forward_mw_real_pole_adjoint(np.ndarray[ double complex, ndim=1, mode="c"] f_lm not None, int L):

        cdef ssht_dl_method_t dl_method = SSHT_DL_RISBO;
        f_mw_r = np.empty([L-1,2*L-1,], dtype=np.float_)
        cdef double f_sp
        ssht_adjoint_mw_forward_sov_sym_real_pole(<double*> np.PyArray_DATA(f_mw_r), &f_sp, <const double complex*> np.PyArray_DATA(f_lm), L, dl_method, 0);
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

def forward(f, int L, int Spin=0, str Method='MW', bint Reality=False):
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
        f_new = np.empty(sample_shape(L,Method=Method), dtype=np.complex_)
        f_new = f + 1j*np.zeros(sample_shape(L,Method=Method), dtype=np.float_)
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


def inverse(flm, int L, int Spin=0, str Method='MW', bint Reality=False):
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

def inverse_adjoint(f, int L, int Spin=0, str Method='MW', bint Reality=False):
    # Checks

    if Method == 'MW_pole':
        if Reality:
            f, f_sp = f
        else:
            f, f_sp, phi_sp = f

    if f.ndim != 2:
      raise ssht_input_error('f must be 2D numpy array')

    if not(Method == 'MW' or Method == 'MWSS'):
        raise ssht_input_error('Method is not recognised, Methods are: MW, MWSS')

    if Spin != 0 and Reality == True :
        raise ssht_spin_error('Reality set to True and Spin is not 0. However, spin signals must be complex.')

    if f.dtype == np.float_ and Reality == False:
        print 'Real signal given but Reality flag is False. Set Reality = True to improve performance'
        f_new = np.empty(sample_shape(L,Method=Method), dtype=np.complex_)
        f_new = f + 1j*np.zeros(sample_shape(L,Method=Method), dtype=np.float_)
        f = f_new

    if f.dtype == np.complex_ and Reality == True:
        print 'Complex signal given but Reality flag is True. Ignoring complex component'
        f_new = np.real(f)
        f = f_new.copy(order='c')
        
    # do correct transform
    if Method == 'MW':
        if Reality:
            flm = ssht_inverse_mw_real_adjoint(f,L)
        else:
            flm = ssht_inverse_mw_complex_adjoint(f,L,Spin)
            
    # if Method == 'MW_pole':
    #     if Reality:
    #         flm = ssht_inverse_mw_real_pole_adjoint(f,f_sp,L)
    #     else:
    #         flm = ssht_inverse_mw_complex_pole_adjoint(f,f_sp,phi_sp,L,Spin)
            
    if Method == 'MWSS':
        if Reality:
            flm = ssht_inverse_mwss_real_adjoint(f,L)
        else:
            flm = ssht_inverse_mwss_complex_adjoint(f,L,Spin)
            
            
    return flm


def forward_adjoint(flm, int L, int Spin=0, str Method='MW', bint Reality=False):
    if flm.ndim != 1:
      raise ssht_input_error('flm must be 1D numpy array')

    if not(Method == 'MW' or Method == 'MWSS'):
        raise ssht_input_error('Method is not recognised, Methods are: MW, MWSS')
    
    if Spin != 0 and Reality == True :
        raise ssht_spin_error('Reality set to True and Spin is not 0. However, spin signals must be complex.')
    
    # do correct transform
    if Method == 'MW':
        if Reality:
            f = ssht_forward_mw_real_adjoint(flm,L)
        else:
            f = ssht_forward_mw_complex_adjoint(flm,L,Spin)
            

    # if Method == 'MW_pole':
    #     if Reality:
    #         f = ssht_forward_mw_real_pole_adjoint(flm,L)
    #     else:
    #         f = ssht_forward_mw_complex_pole_adjoint(flm,L,Spin)

    if Method == 'MWSS':
         if Reality:
             f = ssht_forward_mwss_real_adjoint(flm,L)
         else:
             f = ssht_forward_mwss_complex_adjoint(flm,L,Spin)
                      
    return f
#----------------------------------------------------------------------------------------------------#

# index to ell em and back function

def isqrt(int n):
    cdef int square = 1, delta = 3
    while square < n:
        square += delta
        delta  += 2
    return (delta/2 -1)

def elm2ind( int el, int m):

  return el * el + el + m


def ind2elm(int ind):

  cdef int ell, em
  el = cy_isqrt(ind)
  em = ind - (el)*(el) - (el);

  return el, em

# index to ell em and back function

cdef int cy_isqrt(int n):
    cdef int square = 1, delta = 3
    while square < n:
        square += delta
        delta  += 2
    return (delta/2 -1)

cdef int cy_elm2ind( int el, int m):

  return el * el + el + m


def theta_to_index(double theta, int L, str Method="MW"):
  cdef int p
  cdef np.ndarray[np.float_t, ndim=1] theta_gl_grid, phi_gl_grid

  if Method == 'GL':
    theta_gl_grid, phi_gl_grid = sample_positions(L,Method="GL")

  if Method == 'MW':
    p = int((theta*(2*L-1)/np.pi-1)/2)  # (2.0*t + 1.0) * SSHT_PI / (2.0*L - 1.0)
  if Method == 'MWSS':
    p = int((theta*(2*L)/np.pi)/2)  # 2.0 * t * SSHT_PI / (2.0 * L)
  if Method == 'DH':
    p = int((theta*(4*L)/np.pi-1)/2)  # (2.0*t + 1.0) * SSHT_PI / (4.0*L)
  if Method == 'GL':
    if theta > theta_gl_grid[L-1]:
      p = L-1
      for k in range(theta_gl_grid.size):
        if theta < theta_gl_grid[k]:
          p = k-1
          break
  return p

def phi_to_index(double phi, int L, str Method="MW"):
  cdef int q
  
  if Method == 'MW':
    q = int(phi*(2*L-1)/(2*np.pi))      # 2.0 * p * SSHT_PI / (2.0*L - 1.0)
  if Method == 'MWSS':
    q = int(phi*(2*L)/(2*np.pi))      # 2.0 * p * SSHT_PI / (2.0*L)
  if Method == 'DH':
    q = int(phi*(2*L-1)/(2*np.pi))      # 2.0 * p * SSHT_PI / (2.0*L - 1.0)
  if Method == 'GL':
    q = int(phi*(2*L-1)/(2*np.pi))      # 2.0 * p * SSHT_PI / (2.0*L - 1.0)
  return q

cdef int cy_theta_to_index(double theta, int L, METHOD_TYPE Method_enum):
  cdef int p
  cdef np.ndarray[np.float_t, ndim=1] theta_gl_grid, phi_gl_grid

  if Method_enum == GL:
    theta_gl_grid, phi_gl_grid = sample_positions(L,Method="GL")

  if Method_enum == MW:
    p = int((theta*(2*L-1)/np.pi-1)/2)  # (2.0*t + 1.0) * SSHT_PI / (2.0*L - 1.0)
  if Method_enum == MWSS:
    p = int((theta*(2*L)/np.pi)/2)  # 2.0 * t * SSHT_PI / (2.0 * L)
  if Method_enum == DH:
    p = int((theta*(4*L)/np.pi-1)/2)  # (2.0*t + 1.0) * SSHT_PI / (4.0*L)
  if Method_enum == GL:
    if theta > theta_gl_grid[L-1]:
      p = L-1
      for k in range(theta_gl_grid.size):
        if theta < theta_gl_grid[k]:
          p = k-1
          break
  return p

cdef int cy_phi_to_index(double phi, int L, METHOD_TYPE Method_enum):
  cdef int q
  
  if Method_enum == MW:
    q = int(phi*(2*L-1)/(2*np.pi))      # 2.0 * p * SSHT_PI / (2.0*L - 1.0)
  if Method_enum == MWSS:
    q = int(phi*(2*L)/(2*np.pi))      # 2.0 * p * SSHT_PI / (2.0*L)
  if Method_enum == DH:
    q = int(phi*(2*L-1)/(2*np.pi))      # 2.0 * p * SSHT_PI / (2.0*L - 1.0)
  if Method_enum == GL:
    q = int(phi*(2*L-1)/(2*np.pi))      # 2.0 * p * SSHT_PI / (2.0*L - 1.0)
  return q


def sample_length(int L, Method = 'MW'):
  if not(Method == 'MW' or Method == 'MW_pole' or Method == 'MWSS' or Method == 'DH' or Method == "GL"):
      raise ssht_input_error('Method is not recognised, Methods are: MW, MW_pole, MWSS, DH and GL')

  cdef int n_theta, n_phi
  n_theta, n_phi = sample_shape(L, Method=Method)

  return n_theta*n_phi
# get the shape of the signal on the sphere for different sampling theorems

def sample_shape(int L, str Method = 'MW'):
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

def sample_positions(int L, str Method = 'MW', bint Grid=False):
  if not(Method == 'MW' or Method == 'MW_pole' or Method == 'MWSS' or Method == 'DH' or Method == "GL"):
    raise ssht_input_error('Method is not recognised, Methods are: MW, MW_pole, MWSS, DH and GL')

  cdef int n_theta, n_phi, i
  cdef np.ndarray[np.float_t, ndim=1] thetas, phis

  n_theta, n_phi = sample_shape(L, Method=Method)
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
    phis_grid, thetas_grid = np.meshgrid(phis, thetas)
    return thetas_grid, phis_grid
  else:
    return thetas, phis


def s2_to_cart(theta, phi):
  if theta.shape != phi.shape:
    raise ssht_input_error('theta and phi must be the same shape')

  x = np.sin(theta) * np.cos(phi)
  y = np.sin(theta) * np.sin(phi)
  z = np.cos(theta)
  
  return (x, y, z)

def cart_to_s2(x, y, z):
  if x.shape != y.shape or y.shape != z.shape:
    raise ssht_input_error('x, y and z must be the same shape')

  phi = np.arctan2(y,x)
  theta = np.arctan2(np.sqrt(x*x + y*y),z)

  with np.errstate(invalid='ignore'):
    phi[phi<0] += 2*np.pi 

  return (theta, phi)


def spherical_to_cart(r, theta, phi):
  if theta.shape != r.shape or theta.shape != phi.shape:
    raise ssht_input_error('r, theta and phi must be the same shape.')
  
  x = r * np.sin(theta) * np.cos(phi)
  y = r * np.sin(theta) * np.sin(phi)
  z = r * np.cos(theta)
  
  return (x, y, z)

def cart_to_spherical(x, y, z):
  if x.shape != y.shape or y.shape != z.shape:
    raise ssht_input_error('x, y and z must be the same shape')

  phi = np.arctan2(y,x)
  theta = np.arctan2(np.sqrt(x*x + y*y),z)
  r = np.sqrt(x*x + y*y + z*z)

  with np.errstate(invalid='ignore'):
    phi[phi<0] += 2*np.pi 

  return (r, theta, phi)

def theta_phi_to_ra_dec(theta, phi, bint Degrees=False):
  dec = (theta - np.pi/2)*(-1)
  ra  = phi - np.pi

  if Degrees:
    dec = dec*180/np.pi
    ra  = ra*180/np.pi

  return dec, ra


def ra_dec_to_theta_phi(ra, dec, bint Degrees=False):
  if Degrees:
    dec = dec*np.pi/180
    ra  = ra*np.pi/180

  theta = np.pi/2 - dec
  phi  = ra + np.pi


  return theta, phi

def make_rotation_matrix(list rot):
  cdef np.ndarray[np.float_t, ndim=2] rot_matix = np.empty((3,3), dtype=np.float_)

  cdef double c1 = np.cos(rot[0])
  cdef double s1 = np.sin(rot[0])
  cdef double c2 = np.cos(rot[1])
  cdef double s2 = np.sin(rot[1])
  cdef double c3 = np.cos(rot[2])
  cdef double s3 = np.sin(rot[2])

  rot_matix[0,0] = c3*c1 - c2*s1*s3
  rot_matix[0,1] = c3*s1 + c2*c1*s3
  rot_matix[0,2] = s3*s2
  rot_matix[1,0] = -s3*c1 - c2*s1*c3
  rot_matix[1,1] = -s3*s1 + c2*c1*c3
  rot_matix[1,2] = c3*s2
  rot_matix[2,0] = s2*s1
  rot_matix[2,1] = -s2*c1
  rot_matix[2,2] = c2

  return rot_matix


def rot_cart(x, y, z, list rot):
  if x.shape != y.shape or y.shape != z.shape:
    raise ssht_input_error('x, y and z must be the same shape')

  cdef np.ndarray[np.float_t, ndim=2] rot_matix
  rot_matix = make_rotation_matrix(rot)

  cdef unsigned long i

  x_p = np.empty(x.shape, dtype=np.float_)
  y_p = np.empty(x.shape, dtype=np.float_)
  z_p = np.empty(x.shape, dtype=np.float_)


  for i in range(x.size):
    x_p.flat[i] = x.flat[i]*rot_matix[0,0] + y.flat[i]*rot_matix[0,1] + z.flat[i]*rot_matix[0,2]
    y_p.flat[i] = x.flat[i]*rot_matix[1,0] + y.flat[i]*rot_matix[1,1] + z.flat[i]*rot_matix[1,2]
    z_p.flat[i] = x.flat[i]*rot_matix[2,0] + y.flat[i]*rot_matix[2,1] + z.flat[i]*rot_matix[2,2]

  return x_p, y_p, z_p


def rot_cart_1d(np.ndarray[np.float_t, ndim=1] x, np.ndarray[np.float_t, ndim=1] y, np.ndarray[np.float_t, ndim=1] z, list rot):
  if x.shape[0] != y.shape[0] or y.shape[0] != z.shape[0]:
    raise ssht_input_error('x, y and z must be the same shape')

  cdef np.ndarray[np.float_t, ndim=1] x_p, y_p, z_p
  cdef np.ndarray[np.float_t, ndim=2] rot_matix
  cdef int i, j

  rot_matix = make_rotation_matrix(rot)

  x_p = np.empty((x.shape[0]), dtype=np.float_)
  y_p = np.empty((x.shape[0]), dtype=np.float_)
  z_p = np.empty((x.shape[0]), dtype=np.float_)


  for i in range(x.shape[0]):
      x_p[i] = x[i]*rot_matix[0,0] + y[i]*rot_matix[0,1] + z[i]*rot_matix[0,2]
      y_p[i] = x[i]*rot_matix[1,0] + y[i]*rot_matix[1,1] + z[i]*rot_matix[1,2]
      z_p[i] = x[i]*rot_matix[2,0] + y[i]*rot_matix[2,1] + z[i]*rot_matix[2,2]

  return x_p, y_p, z_p

def rot_cart_2d(np.ndarray[np.float_t, ndim=2] x, np.ndarray[np.float_t, ndim=2] y, np.ndarray[np.float_t, ndim=2] z, list rot):
  if x.shape[0] != y.shape[0] or y.shape[0] != z.shape[0] or x.shape[1] != y.shape[1] or y.shape[1] != z.shape[1]:
    raise ssht_input_error('x, y and z must be the same shape')

  cdef np.ndarray[np.float_t, ndim=2] x_p, y_p, z_p, rot_matix
  cdef int i, j

  rot_matix = make_rotation_matrix(rot)

  x_p = np.empty((x.shape[0], x.shape[1]), dtype=np.float_)
  y_p = np.empty((x.shape[0], x.shape[1]), dtype=np.float_)
  z_p = np.empty((x.shape[0], x.shape[1]), dtype=np.float_)


  for i in range(x.shape[0]):
    for j in range(x.shape[1]):
      x_p[i,j] = x[i,j]*rot_matix[0,0] + y[i,j]*rot_matix[0,1] + z[i,j]*rot_matix[0,2]
      y_p[i,j] = x[i,j]*rot_matix[1,0] + y[i,j]*rot_matix[1,1] + z[i,j]*rot_matix[1,2]
      z_p[i,j] = x[i,j]*rot_matix[2,0] + y[i,j]*rot_matix[2,1] + z[i,j]*rot_matix[2,2]

  return x_p, y_p, z_p


#----------------------------------------------------------------------------------------------------#

# Plotting functions

def plot_sphere(f, int L, str Method='MW', bint Close=True, bint Parametric=False, list Parametric_Saling=[0.0,0.5], \
                     Output_File=None, bint Show=True, bint Color_Bar=True, Units=None, Color_Range=None, \
                     int Axis=True): 

    # add ability to choose color bar min max
    # and sort out shapes of the plots
    
    if Method == 'MW_pole':
        if len(f) == 2:
            f, f_sp = f
        else:
            f, f_sp, phi_sp = f

    (thetas, phis) = sample_positions(L, Method=Method, Grid=True);

    if (thetas.size != f.size):
        ssht_input_error('Band limit L deos not match that of f')

    f_plot = f.copy()

    f_max = f_plot.max()
    f_min = f_plot.min()

    if Color_Range is None:
        f_max = f_plot.max()
        f_min = f_plot.min()
        vmin = f_min
        vmax = f_max
    else:
        vmin = Color_Range[0]
        vmax = Color_Range[1]
        f_plot[f_plot<Color_Range[0]] = Color_Range[0]
        f_plot[f_plot>Color_Range[1]] = Color_Range[1]
        f_plot[f_plot == -1.56E30] = np.nan

    #% Compute position scaling for parametric plot.
    if Parametric:
        f_normalised = (f_plot - vmin/(vmax - vmin))*Parametric_Saling[1]+Parametric_Saling[0]


    # % Close plot.
    if Close:
        (n_theta,n_phi) = sample_shape(L,Method=Method)
        f_plot = np.insert(f_plot,n_phi,f[:,0],axis=1)
        if Parametric:
          f_normalised = np.insert(f_normalised,n_phi,f_normalised[:,0],axis=1)
        thetas = np.insert(thetas,n_phi,thetas[:,0],axis=1)
        phis = np.insert(phis,n_phi,phis[:,0],axis=1)



    
    # % Compute location of vertices.
    if Parametric:
        (x, y, z) = spherical_to_cart(f_normalised, thetas, phis)
    else:
        (x, y, z) = s2_to_cart(thetas, phis)
        

    # % Plot.
    fig = plt.figure(figsize=plt.figaspect(1.1))
    gs = gridspec.GridSpec(2,1, height_ratios=[10,0.5])
    ax = fig.add_subplot(gs[0], projection='3d')
    ax_cbar = fig.add_subplot(gs[1])
    norm = colors.Normalize()
    surf = ax.plot_surface(x, y, z,  rstride=1, cstride=1, facecolors=cm.jet(norm(f_plot)))
    if not Axis:
        ax.set_axis_off()


    if Color_Bar:
        cmap = cm.jet
        norm = colors.Normalize(vmin=vmin, vmax=vmax)
        cb1 = colorbar.ColorbarBase(ax_cbar, cmap=cmap,
                                    norm=norm,
                                    orientation='horizontal')
        if Units != None:
            cb1.set_label(Units)

   
    # output (to file and screan) 
    if Output_File != None:
        plt.savefig(Output_File)
    if Show:
        plt.show()
        
    return



def mollweide_projection_work(np.ndarray[ double, ndim=2, mode="c"] f, int L, int resolution=500, rot=None,\
                        zoom_region=[np.sqrt(2.0)*2,np.sqrt(2.0)], METHOD_TYPE Method_enum=MW):

  if not(len(zoom_region)==2 or len(zoom_region)==4):
    raise ssht_input_error('zoom_region must be a python list of length 2 or 4')

  if rot is None or len(rot) > 1:
    dummy = 0.0
  else:
    dummy = <double> rot[0]

  cdef double rot_angle = dummy

  cdef int  Nx, Ny 
  cdef int p, q, i, j

  cdef np.ndarray[np.float_t, ndim=1] x1, y1
  cdef np.ndarray[np.float_t, ndim=2] x, y, aux_angle, dec, ra, theta, phi
  cdef np.ndarray[np.float_t, ndim=2] xx, yy, zz, xx_p, yy_p, zz_p, f_plot, mask



  if len(zoom_region) == 2:
    Nx = resolution*zoom_region[0]/zoom_region[1]
    Ny = resolution
    x1 = np.linspace(-zoom_region[0],zoom_region[0],num=Nx,dtype=np.float_)
    y1 = np.linspace(-zoom_region[1],zoom_region[1],num=Ny,dtype=np.float_)
  if len(zoom_region) == 4:
    Nx = resolution*(zoom_region[1]-zoom_region[0])/(zoom_region[3]-zoom_region[2])
    Ny = resolution
    x1 = np.linspace(zoom_region[0],zoom_region[1],num=Nx,dtype=np.float_)
    y1 = np.linspace(zoom_region[2],zoom_region[3],num=Ny,dtype=np.float_)

  x, y = np.meshgrid(x1,y1)


  aux_angle = np.arcsin(y/np.sqrt(2.0))
  dec = np.arcsin((2*aux_angle + np.sin(2*aux_angle))/np.pi)
  ra = np.pi*x/(2*np.sqrt(2.0)*np.cos(aux_angle))

  theta, phi = ra_dec_to_theta_phi(ra, dec)

  for i from 0 <= i < Ny:
    for j from 0 <= j < Nx:
      if phi[i,j]<0 or phi[i,j]>2*np.pi:
        phi[i,j] = np.nan

  if not(rot is None):
    if len(rot) == 1:
      for i from 0 <= i < Ny:
        for j from 0 <= j < Nx:
          if not(np.isnan(phi[i,j])):
            phi[i,j] += rot_angle
            if phi[i,j] > 2*np.pi:
              phi[i,j] -= 2*np.pi
            if phi[i,j] < 0:
              phi[i,j] += 2*np.pi

    elif len(rot) == 3:
      xx, yy, zz = s2_to_cart(theta, phi)
      xx_p, yy_p, zz_p = rot_cart_2d(xx, yy, zz, rot)
      theta, phi = cart_to_s2(xx_p, yy_p, zz_p)

  f_plot = np.empty([Ny,Nx],dtype=np.float_)
  mask   = np.empty([Ny,Nx],dtype=np.float_)

  mask.fill(np.nan)

  # bin the values
  for i from 0 <= i < Ny:
    for j from 0 <= j < Nx:
      if np.isnan(theta[i,j]) or np.isnan(phi[i,j]):
        f_plot[i,j] = np.nan
      else:
        p = cy_phi_to_index(theta[i,j], L, Method_enum)
        q = cy_phi_to_index(phi[i,j], L, Method_enum)


        if np.isnan(f[p,q]):
          f_plot[i,j] = np.nan
          mask[i,j] = 0.0
        else:
          f_plot[i,j] = f[p,q]

  return f_plot, mask


def mollweide_projection(f, int L, int resolution=500, rot=None,\
                        zoom_region=[np.sqrt(2.0)*2,np.sqrt(2.0)], str Method="MW"):

  cdef int i, j, n_phi, n_theta
  cdef np.ndarray[np.float_t, ndim=2] f_real, f_imag
  cdef METHOD_TYPE Method_enum

  if Method=="MW":
    Method_enum=MW
  elif Method=="MWSS":
    Method_enum = MWSS
  elif Method=="DH":
    Method_enum = DH
  elif Method=="GL":
    Method_enum = GL
  else:
    raise ssht_input_error('Method is not recognised, Methods are: MW, MWSS, DH and GL')

  if not isinstance(f, np.ndarray):
    raise TypeError("Input not a ndarray")

  if f.ndim != 2:
    raise ssht_input_error("f must have ndim = 2")


  if f.dtype == np.float64:
    return mollweide_projection_work(f, L, resolution=resolution, rot=rot,\
                        zoom_region=zoom_region, Method_enum=Method_enum)
  elif f.dtype == complex:
    n_theta, n_phi = sample_shape(L,Method=Method)
    f_real = np.empty((n_theta,n_phi), dtype=np.float_)
    f_imag = np.empty((n_theta,n_phi), dtype=np.float_)
    for i in range(n_theta):
      for j in range(n_phi):
        f_real[i,j] = f[i,j].real
        f_imag[i,j] = f[i,j].imag

    f_real_plot, mask_real = mollweide_projection_work(f_real, L, resolution=resolution, rot=rot,\
                        zoom_region=zoom_region, Method_enum=Method_enum)
    f_imag_plot, mask_imag = mollweide_projection_work(f_imag, L, resolution=resolution, rot=rot,\
                        zoom_region=zoom_region, Method_enum=Method_enum)
    return f_real_plot, mask_real, f_imag_plot, mask_imag
  else:
    raise ssht_input_error("f dtype must be float or complex")

def orthographic_projection_work(np.ndarray[ double, ndim=2, mode="c"] f, int L, int resolution=500, rot=None,\
                        float zoom_region=np.pi/2,  METHOD_TYPE Method_enum=MW):

  if rot is None or len(rot) > 1:
    dummy = 0.0
  else:
    dummy = <double> rot[0]

  if Method_enum==MW:
    Method="MW"
  elif Method_enum==MWSS:
    Method_enum = "MWSS"
  elif Method_enum==DH:
    Method = "DH"
  elif Method_enum==GL:
    Method = "GL"

  cdef double rot_angle = dummy

  cdef np.ndarray[np.float_t, ndim=2] theta, phi
  cdef np.ndarray[np.float_t, ndim=2] x_project, y_project
  cdef np.ndarray[np.float_t, ndim=2] xx, yy, zz, xx_p, yy_p, zz_p
  cdef np.ndarray[np.float_t, ndim=2] ortho_proj_north, mask_north, ortho_proj_south, mask_south
  cdef np.ndarray[np.float_t, ndim=2] rot_matix

  cdef list rot_list

  cdef np.ndarray[np.int_t,   ndim=2] n_points_north, n_points_south

  cdef float x_pos, y_pos, z_pos, x_p_pos, y_p_pos, z_p_pos, rho, theta_pos, phi_pos, half_box_len
  
  cdef int n_theta, n_phi, n_theta_north, n_theta_south, i, j, p, q, i_rot, j_rot

  half_box_len = np.sin(zoom_region)

  n_theta, n_phi = sample_shape(L, Method=Method)
  n_theta_north = n_theta/2
  n_theta_south = n_theta-n_theta_north

  theta, phi = sample_positions(L, Grid=True, Method=Method)

  if not(rot is None):
    if len(rot) == 1:
      for i from 0 <= i < n_theta:
        for j from 0 <= j < n_phi:
          phi[i,j] += rot_angle
          if phi[i,j] > 2*np.pi:
            phi[i,j] -= 2*np.pi
          if phi[i,j] < 0:
            phi[i,j] += 2*np.pi

    elif len(rot) == 3:
      xx, yy, zz = s2_to_cart(theta, phi)
      xx_p, yy_p, zz_p = rot_cart_2d(xx, yy, zz, rot)
      theta, phi = cart_to_s2(xx_p, yy_p, zz_p)

  # do projection
  x_project = -np.sin(theta)*np.cos(phi)
  y_project = -np.sin(theta)*np.sin(phi)

  ortho_proj_north = np.zeros((resolution, resolution), dtype=float)
  ortho_proj_south = np.zeros((resolution, resolution), dtype=float)

  n_points_north = np.zeros((resolution, resolution), dtype=int)
  n_points_south = np.zeros((resolution, resolution), dtype=int)

  mask_north   = np.empty([resolution,resolution],dtype=np.float_)
  mask_south   = np.empty([resolution,resolution],dtype=np.float_)

  mask_north.fill(np.nan)
  mask_south.fill(np.nan)

  for i in range(n_theta):
    for j in range(n_phi):

      if x_project[i,j] < half_box_len and x_project[i,j] > -half_box_len and y_project[i,j] < half_box_len and y_project[i,j] > -half_box_len:
        p = int(resolution*(x_project[i,j]+half_box_len)/(2.0*half_box_len))
        q = int(resolution*(y_project[i,j]+half_box_len)/(2.0*half_box_len))

        if theta[i,j] < np.pi/2:        
          if np.isnan(f[i,j]):
            ortho_proj_north[p,q] = np.nan
            mask_north[p,q] = 0.0
          else:
            ortho_proj_north[p,q] += f[i,j]
            n_points_north[p,q]   += 1
        else:
          if np.isnan(f[i,j]):
            ortho_proj_south[p,q] = np.nan
            mask_south[p,q] = 0.0
          else:
            ortho_proj_south[p,q] += f[i,j]
            n_points_south[p,q]   += 1

  if not(rot is None):
    if len(rot) == 1:
      rot_list = [rot[0],0.0,0.0]
    elif len(rot) == 3:
      rot_list = [-rot[2],-rot[1],-rot[0]]
    rot_matix = make_rotation_matrix(rot_list)

  for i in range(resolution):
    for j in range(resolution):
      if n_points_north[i,j] == 0:
        if orthographic_projection_index_to_length(i,j,resolution, half_box_len) < 1.0:
          x_pos     = -(2.*half_box_len*(<float>i)/<float>resolution -half_box_len)
          y_pos     = -(2.*half_box_len*(<float>j)/<float>resolution -half_box_len)
          rho       = np.sqrt(x_pos*x_pos + y_pos*y_pos) 
          theta_pos = np.arcsin(rho)
          phi_pos   = np.arctan2(y_pos,x_pos)

          # perform rotation
          if rot is not None:
            z_pos   = np.cos(theta_pos)
            x_p_pos = x_pos*rot_matix[0,0] + y_pos*rot_matix[0,1] + z_pos*rot_matix[0,2]
            y_p_pos = x_pos*rot_matix[1,0] + y_pos*rot_matix[1,1] + z_pos*rot_matix[1,2]
            z_p_pos = x_pos*rot_matix[2,0] + y_pos*rot_matix[2,1] + z_pos*rot_matix[2,2]

            theta_pos = np.arctan2(np.sqrt(x_p_pos*x_p_pos + y_p_pos*y_p_pos),z_p_pos)
            phi_pos   = np.arctan2(y_p_pos,x_p_pos)

          if phi_pos < 0:
            phi_pos += 2*np.pi
          if phi_pos > 2*np.pi:
            phi_pos -= 2*np.pi
          
          p = cy_phi_to_index(theta_pos, L, Method_enum)
          q = cy_phi_to_index(phi_pos, L, Method_enum)
          if np.isnan(f[p,q]):
            ortho_proj_north[i,j] = np.nan
            mask_north[i,j] = 0.0
          else:
            ortho_proj_north[i,j] = f[p,q]
            n_points_north[i,j]   += 1
        else:
          ortho_proj_north[i,j] = np.nan
      else:
        ortho_proj_north[i,j] = ortho_proj_north[i,j]/<float>n_points_north[i,j]

      if n_points_south[i,j] == 0:
        if orthographic_projection_index_to_length(i,j,resolution, half_box_len) < 1.0:
          x_pos     = -(2.*half_box_len*(<float>i)/<float>resolution -half_box_len)
          y_pos     = -(2.*half_box_len*(<float>j)/<float>resolution -half_box_len)
          rho       = np.sqrt(x_pos*x_pos + y_pos*y_pos) 
          theta_pos = np.pi-np.arcsin(rho)
          phi_pos   = np.arctan2(y_pos,x_pos)

          # perform rotation
          if rot is not None:
            z_pos   = np.cos(theta_pos)
            x_p_pos = x_pos*rot_matix[0,0] + y_pos*rot_matix[0,1] + z_pos*rot_matix[0,2]
            y_p_pos = x_pos*rot_matix[1,0] + y_pos*rot_matix[1,1] + z_pos*rot_matix[1,2]
            z_p_pos = x_pos*rot_matix[2,0] + y_pos*rot_matix[2,1] + z_pos*rot_matix[2,2]

            theta_pos = np.arctan2(np.sqrt(x_p_pos*x_p_pos + y_p_pos*y_p_pos),z_p_pos)
            phi_pos   = np.arctan2(y_p_pos,x_p_pos)

          if phi_pos < 0:
            phi_pos += 2*np.pi
          if phi_pos > 2*np.pi:
            phi_pos -= 2*np.pi
          
          p = cy_phi_to_index(theta_pos, L, Method_enum)
          q = cy_phi_to_index(phi_pos, L, Method_enum)
          if np.isnan(f[p,q]):
            ortho_proj_south[i,j] = np.nan
            mask_south[i,j] = 0.0
          else:
            ortho_proj_south[i,j] = f[p,q]
            n_points_south[i,j]   += 1
        else:
          ortho_proj_south[i,j] = np.nan
      else:
        ortho_proj_south[i,j] = ortho_proj_south[i,j]/<float>n_points_south[i,j]

  return ortho_proj_north, mask_north, ortho_proj_south, mask_south


cdef float orthographic_projection_index_to_length(int i, int j, int N, float half_box_len):
  return (2.*half_box_len*<float>i/<float>N -half_box_len)*(2.*half_box_len*<float>i/<float>N -half_box_len) \
          + (2.*half_box_len*<float>j/<float>N -half_box_len)*(2.*half_box_len*<float>j/<float>N -half_box_len)

def orthographic_projection(f, int L, int resolution=500, rot=None,\
                        float zoom_region=np.pi/2, str Method="MW"):

  cdef int i, j, n_phi, n_theta
  cdef np.ndarray[np.float_t, ndim=2] f_real, f_imag
  cdef METHOD_TYPE Method_enum

  if Method=="MW":
    Method_enum=MW
  elif Method=="MWSS":
    Method_enum = MWSS
  elif Method=="DH":
    Method_enum = DH
  elif Method=="GL":
    Method_enum = GL
  else:
    raise ssht_input_error('Method is not recognised, Methods are: MW, MWSS, DH and GL')

  if not isinstance(f, np.ndarray):
    raise TypeError("Input not a ndarray")

  if f.ndim != 2:
    raise ssht_input_error("f must have ndim = 2")


  if f.dtype == np.float64:
    return orthographic_projection_work(f, L, resolution=resolution, rot=rot,\
                        zoom_region=zoom_region, Method_enum=Method_enum)
  elif f.dtype == complex:
    n_theta, n_phi = sample_shape(L,Method=Method)
    f_real = np.empty((n_theta,n_phi), dtype=np.float_)
    f_imag = np.empty((n_theta,n_phi), dtype=np.float_)
    for i in range(n_theta):
      for j in range(n_phi):
        f_real[i,j] = f[i,j].real
        f_imag[i,j] = f[i,j].imag

    north_plot_real, mask_north_real, south_plot_real, mask_south_real\
                         = orthographic_projection_work(f_real, L, resolution=resolution, rot=rot,\
                        zoom_region=zoom_region, Method_enum=Method_enum)
    north_plot_imag, mask_north_imag, south_plot_imag, mask_south_imag\
                        = orthographic_projection_work(f_imag, L, resolution=resolution, rot=rot,\
                        zoom_region=zoom_region, Method_enum=Method_enum)
    return north_plot_real, mask_north_real, south_plot_real, mask_south_real,\
           north_plot_imag, mask_north_imag, south_plot_imag, mask_south_imag
  else:
    raise ssht_input_error("f dtype must be float or complex")



def mollweide_coords_s2_to_xy(thetas, phis):
  # % ssht_mollweide - Compute Mollweide projection
  # %
  # % Compute Mollweide projection of spherical coordinates.
  # %
  # % Usage is given by
  # %
  # %   (x,y) = mollweide_coords(thetas, phis)
  # %
  # % where thetas and phis are spherical coordinates and x and y are the
  # % projected Mollweide coordinates.

  MAX_ITERATIONS = 1e5;
  TOL = 1e-10;

  #% Convert theta to longitude.
  thetas, phis = theta_phi_to_ra_dec(thetas,phis)

  t = thetas
  cdef bint inaccurate = True
  cdef int count = 0

  while(inaccurate):
    count += 1
    dt = (t + np.sin(t) - np.pi*np.sin(thetas)) / (1 + np.cos(t));
    t = t - dt;
    if(np.max(np.abs(dt)) < TOL or count > MAX_ITERATIONS):
      inaccurate = False

  t = t/2;
  x = 2 * np.sqrt(2) / np.pi * phis * np.cos(t);
  y = np.sqrt(2) * np.sin(t);

  return (x, y)

#----------------------------------------------------------------------------------------------------#

# dl functions

def dl_beta_recurse(np.ndarray[ double, ndim=2, mode="c"] dl not None, double beta, int L,\
            int el, np.ndarray[ double, ndim=1, mode="c"] sqrt_tbl not None,\
            np.ndarray[ double, ndim=1, mode="c"] signs not None):

  cdef ssht_dl_size_t dl_size=SSHT_DL_HALF
  
  ssht_dl_beta_risbo_half_table(<double*> np.PyArray_DATA(dl), beta, L, dl_size,\
                                           el, <double*> np.PyArray_DATA(sqrt_tbl), <double*> np.PyArray_DATA(signs));

  return dl

def generate_dl(double beta, int L):

  cdef np.ndarray[np.float_t, ndim=3] dl_array = np.empty((L, 2*L-1, 2*L-1), dtype=np.float_)
  cdef np.ndarray[np.float_t, ndim=2] dl_dummy = np.zeros((2*L-1, 2*L-1), dtype=np.float_)

  cdef np.ndarray[np.float_t, ndim=1] sqrt_tbl = np.sqrt(np.arange(0,2*(L-1)+1, dtype=np.float_))
  cdef np.ndarray[np.float_t, ndim=2] signs = np.ones((L+1,1), dtype=np.float_)

  cdef int i, el

  for i in range(1,L+1,2):
    signs[i] = -1 

  
  cdef ssht_dl_size_t dl_size=SSHT_DL_HALF
  
  # do recursion
  #el = 0 first
  ssht_dl_beta_risbo_half_table(<double*> np.PyArray_DATA(dl_dummy), beta, L, dl_size,\
                                           0, <double*> np.PyArray_DATA(sqrt_tbl), <double*> np.PyArray_DATA(signs));
  dl_array[0,:,:] = dl_dummy

  for el in range(1,L):
    ssht_dl_beta_risbo_half_table(<double*> np.PyArray_DATA(dl_dummy), beta, L, dl_size,\
                                           el, <double*> np.PyArray_DATA(sqrt_tbl), <double*> np.PyArray_DATA(signs));
    dl_array[el,:,:] = dl_dummy

  return dl_array


def generate_exp_array(double x, int L):

  cdef np.ndarray[complex, ndim=1] exp_array = np.empty((2*L-1), dtype=np.complex_)

  cdef int i = 0, m

  for m in range(-L+1,L):
    exp_array[i] = np.exp(-1j*<double>m*x)
    i += 1

  return exp_array

#----------------------------------------------------------------------------------------------------#

# rotation functions


def rotate_flms(np.ndarray[ double complex, ndim=1, mode="c"] f_lm not None,\
                      double alpha, double beta, double gamma, int L, dl_array_in=None,\
                      M_in=None, bint Axisymmetric=False, bint Keep_dl=False):

  cdef np.ndarray[np.float_t, ndim=3] dl_array
  cdef np.ndarray[complex, ndim=1] alpha_array, gamma_array 
  cdef np.ndarray[complex, ndim=1] f_lm_rotated
  cdef int M, el, m

  cdef int index = 0, ind, n_max, n
  cdef complex Dlmn=0

  if dl_array_in == None:
    dl_array = generate_dl(beta, L)
  else:
    dl_array = dl_array_in

  alpha_array = generate_exp_array(alpha, L)
  gamma_array = generate_exp_array(gamma, L)

  if M_in==None:
    M = L
  else:
    M = M_in



  f_lm_rotated = np.zeros((L*L), dtype=complex)

  for el in range(L):
    for m in range(-el,el+1):
        if Axisymmetric:
            n_max = 0
        else:
            n_max = min(el, M-1)

        for n in range(-n_max,n_max+1):
            Dlmn =  <complex> alpha_array[m+L-1] * <complex> dl_array[el,m+L-1,n+L-1]\
                   * <complex> gamma_array[n+L-1] # not sure about
            if Axisymmetric:
                ind = el;
            else:
                ind = cy_elm2ind(el,n);

            f_lm_rotated[index] = <double complex> f_lm_rotated[index] + \
                <complex> Dlmn * <complex> f_lm[ind];

        index = index + 1;

  if Keep_dl:
    return f_lm_rotated, dl_array
  else:
    return f_lm_rotated


#----------------------------------------------------------------------------------------------------#

# other usefull functions

def guassian_smoothing(np.ndarray[ double complex, ndim=1, mode="c"] f_lm not None, int L, sigma_in=None, bl_in = None):

  cdef double sigma
  cdef np.ndarray[ double, ndim=1, mode="c"] bl
  cdef np.ndarray[ complex, ndim=1, mode="c"] fs_lm = np.empty([L*L], dtype=complex)
  cdef int index, el, m


  if sigma_in is None and bl_in is None:
    raise ssht_input_error('One of sigma or bl must be set')

  if sigma_in is not None and bl_in is not None:
    raise ssht_input_error('Only one of sigma or bl can be set')


  if sigma_in is not None:
    sigma = sigma_in
    bl = np.empty([L,], dtype=float)
    for el in range(L):
      bl[el] = np.exp(-<float>el*<float>el*sigma*sigma)

  if bl_in is not None:
    bl = bl_in

  for el in range(L):
    for m in range(-el,el+1):
      index = cy_elm2ind(el, m)
      fs_lm[index] = f_lm[index]*bl[el]

  return fs_lm