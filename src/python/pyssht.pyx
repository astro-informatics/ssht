# import both numpy and the Cython declarations for numpy
import numpy as np
cimport numpy as np

import matplotlib.pyplot as plt
from matplotlib import cm, colors, colorbar, gridspec
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.axes_grid1 import make_axes_locatable

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

def forward(f, int L, Spin=0, Method='MW', Reality=False):
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


def inverse(flm, L, Spin=0, Method='MW', Reality=False):
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

def inverse_adjoint(f, int L, Spin=0, Method='MW', Reality=False):
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


def forward_adjoint(flm, L, Spin=0, Method='MW', Reality=False):
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
  el = isqrt(ind)
  em = ind - (el)*(el) - (el);

  return el, em

def sample_length(int L, Method = 'MW'):
  if not(Method == 'MW' or Method == 'MW_pole' or Method == 'MWSS' or Method == 'DH' or Method == "GL"):
      raise ssht_input_error('Method is not recognised, Methods are: MW, MW_pole, MWSS, DH and GL')

  n_theta, n_phi = sample_shape(L, Method=Method)

  return n_theta*n_phi
# get the shape of the signal on the sphere for different sampling theorems

def sample_shape(int L, Method = 'MW'):
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

def sample_positions(int L, Method = 'MW', Grid=False):
  if not(Method == 'MW' or Method == 'MW_pole' or Method == 'MWSS' or Method == 'DH' or Method == "GL"):
    raise ssht_input_error('Method is not recognised, Methods are: MW, MW_pole, MWSS, DH and GL')

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
    phis, thetas = np.meshgrid(phis, thetas)

  return thetas, phis


def s2_to_cart(theta, phi):
  if theta.shape != phi.shape:
    raise ssht_input_error('theta and phi must be the same shape')

  x = np.sin(theta) * np.cos(phi)
  y = np.sin(theta) * np.sin(phi)
  z = np.cos(theta)
  return (x, y, z)

def spherical_to_cart(r, theta, phi):
  if theta.shape != r.shape or theta.shape != phi.shape:
    raise ssht_input_error('r, theta and phi must be the same shape.')
  
  x = r * np.sin(theta) * np.cos(phi)
  y = r * np.sin(theta) * np.sin(phi)
  z = r * np.cos(theta)
  
  return (x, y, z)


def theta_phi_to_ra_dec(theta, phi, Degrees=False):
  dec = (theta - np.pi/2)*(-1)
  ra  = phi - np.pi

  if Degrees:
    dec = dec*180/np.pi
    ra  = ra*180/np.pi

  return dec, ra
#----------------------------------------------------------------------------------------------------#

# Plotting functions

def plot_sphere(f, L, Method='MW', Close=True, Parametric=False, Parametric_Saling=[0.0,0.5], \
                     Output_File=None, Show=True,  Color_Bar=True, Units=None, Color_Range=None, \
                     Axis=True): # add int L

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

def plot_mollweide(f, L, Method="MW", Close=True):

  theta, phi = sample_positions(L, Method=Method, Grid=True)
  x, y = mollweide_coords_s2_to_xy(theta, phi)

  f_plot = f.copy()
  (n_theta,n_phi) = sample_shape(L,Method=Method)


  if Close:
        f_plot = np.insert(f_plot, n_phi, f_plot[:,0], axis=1)
        x      = np.insert(x,      n_phi, x[:,0],      axis=1)
        y      = np.insert(y,      n_phi, y[:,0],      axis=1)



  N_max = 2*128*128
  if f_plot.size > N_max:
    step = f_plot.size/N_max
    
    x_av = x[range(0,x.size,step),range(0,n_phi,step)]
    y_av = y[range(0,n_theta,step),range(0,n_phi,step)]
    f_plot_av = f_plot[range(0,n_theta,step),range(0,n_phi,step)]
    for i in range(1,step):
      x_av += x[range(i,n_theta,step),range(i,n_phi,step)]
      y_av += y[range(i,n_theta,step),range(i,n_phi,step)]
      f_plot_av += f_plot[range(i,n_theta,step),range(i,n_phi,step)]
    x = x_av/step
    y = y_av/step
    f_plot = f_plot_av/step


  x = x.flatten()
  y = y.flatten()
  f_plot = f_plot.flatten()


  levels = np.arange(f_plot.min(), f_plot.max(), (f_plot.max()-f_plot.min())/50)

  m =  plt.tricontourf(x, y, f_plot, levels, cmap=cm.jet)

  return m

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

  dl_array = np.empty((L, 2*L-1, 2*L-1), dtype=np.float_)
  dl_dummy = np.zeros((2*L-1, 2*L-1), dtype=np.float_)

  sqrt_tbl = np.sqrt(np.arange(0,2*(L-1)+1, dtype=np.float_))
  signs = np.ones((L+1,1), dtype=np.float_)
  for i in range(1,L+1,2):
    signs[i] = -1 


  cdef int stride = np.PyArray_STRIDE(dl_array, 2)  
  strides = np.PyArray_STRIDES(dl_array)  

  
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


def generate_exp_array(double x, L):

  exp_array = np.empty((2*L-1), dtype=np.complex_)

  cdef int i = 0
  for m in range(-L+1,L):
    exp_array[i] = np.exp(-1j*m*x)
    i += 1

  return exp_array

#----------------------------------------------------------------------------------------------------#

# rotation functions


def rotate_flms(np.ndarray[ double complex, ndim=1, mode="c"] f_lm not None,\
                      double alpha, double beta, double gamma, int L, dl_array=None,\
                      M=None, Axisymmetric=False, Keep_dl=False):

  if dl_array == None:
    dl_array = generate_dl(beta, L)

  alpha_array = generate_exp_array(alpha, L)
  gamma_array = generate_exp_array(gamma, L)

  if M==None:
    M = L

  cdef int index = 0, ind;
  cdef double complex Dlmn=0

  f_lm_rotated = np.zeros((L*L), dtype=np.complex_)

  for el in range(L):
    for m in range(-el,el+1):
        if Axisymmetric:
            n_max = 0
        else:
            n_max = min(el, M-1)

        for n in range(-n_max,n_max+1):
            Dlmn =  <double complex> alpha_array[m+L-1] * <double complex> dl_array[el,m+L-1,n+L-1]\
                   * <double complex> gamma_array[n+L-1] # not sure about
            if Axisymmetric:
                ind = el;
            else:
                ind = elm2ind(el,n);

            f_lm_rotated[index] = <double complex> f_lm_rotated[index] + \
                <double complex> Dlmn * <double complex> f_lm[ind];

        index = index + 1;

  if Keep_dl:
    return f_lm_rotated, dl_array
  else:
    return f_lm_rotated

