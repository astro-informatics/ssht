# cython: language_level=3

# import both numpy and the Cython declarations for numpy
import numpy as np
cimport numpy as np
cimport cython

from libc.math cimport log, exp, sqrt, atan2, cos, sin, asin, atan, tan
from scipy.special import factorial
from scipy.interpolate import interp2d
from pyssht.exceptions import ssht_input_error

cdef double pi=3.1415926535897932384626433832795028841971

#----------------------------------------------------------------------------------------------------#

cdef enum METHOD_TYPE:
  MW, MWSS, DH, GL

cdef enum POLAR_PROJECTION_TYPE:
  GP, OP, SP

cdef enum EQUATORIAL_PROJECTION_TYPE:
  SINE, MERCATOR


#----------------------------------------------------------------------------------------------------#
cdef extern from "ssht/ssht.h":

        double ssht_sampling_mw_t2theta(int t, int L) #  adjoints
        double ssht_sampling_mw_p2phi(int p, int L)
        int ssht_sampling_mw_n(int L)
        int ssht_sampling_mw_ntheta(int L)
        int ssht_sampling_mw_nphi(int L)
        double* ssht_dl_calloc(int L, ssht_dl_size_t dl_size)
        ctypedef enum ssht_dl_size_t:
                SSHT_DL_QUARTER_EXTENDED, SSHT_DL_HALF, SSHT_DL_FULL

        void ssht_dl_beta_risbo_full_table(double *dl, double beta, int L,
                                           ssht_dl_size_t dl_size,
                                           int el, double *sqrt_tbl)

        void ssht_dl_beta_risbo_half_table(double *dl, double beta, int L,
                                           ssht_dl_size_t dl_size,
                                           int el, double *sqrt_tbl, double *signs)
        void ssht_dl_beta_kostelec_halfline_table(double *dlm1p1_line,
                                                  double *dl_line,
					                                        double beta, int L, int mm,
                                                  int el, double *sqrt_tbl,
                                                  double *signs)

# I included
        ctypedef enum  ssht_dl_method_t:
                SSHT_DL_RISBO, SSHT_DL_TRAPANI
        void ssht_core_mw_forward_sov_conv_sym(double complex *flm, const double complex *f,
				       int L, int spin,
				       ssht_dl_method_t dl_method,
				       int verbosity)
        void ssht_core_mw_inverse_sov_sym(double complex *f, const double complex *flm,
                                  int L, int spin,
                                  ssht_dl_method_t dl_method,
                                  int verbosity)
        void ssht_core_mw_inverse_sov_sym_real(double *f, const double complex *flm,
                                  int L,
                                  ssht_dl_method_t dl_method,
                                  int verbosity)
        void ssht_core_mw_forward_sov_conv_sym_real(double complex *flm, const double *f,
                                  int L,
                                  ssht_dl_method_t dl_method,
                                  int verbosity)
        void ssht_core_mw_inverse_sov_sym_pole(double complex *f,
                                  double complex *f_sp, double *phi_sp,
                                  const double complex *flm,
                                  int L, int spin,
                                  ssht_dl_method_t dl_method,
                                  int verbosity)
        void ssht_core_mw_inverse_sov_sym_real_pole(double *f,
                                  double *f_sp,
                                  const double complex *flm,
                                  int L,
                                  ssht_dl_method_t dl_method,
                                  int verbosity)
        void ssht_core_mw_forward_sov_conv_sym_pole(double complex *flm, const double complex *f, # change to ss 
                                  double complex f_sp, double phi_sp,
                                  int L, int spin,
                                  ssht_dl_method_t dl_method,
                                  int verbosity)
        void ssht_core_mw_forward_sov_conv_sym_real_pole(double complex *flm,
                                  const double *f,
                                  double f_sp,
                                  int L,
                                  ssht_dl_method_t dl_method,
                                  int verbosity)
        void ssht_core_mw_inverse_sov_sym_ss(double complex *f, const double complex *flm,
              int L, int spin,
              ssht_dl_method_t dl_method,
              int verbosity)
        void ssht_core_mw_inverse_sov_sym_ss_real(double *f, const double complex *flm,
              int L,
              ssht_dl_method_t dl_method,
              int verbosity)
        void ssht_core_mw_forward_sov_conv_sym_ss(double complex *flm, const double complex *f,
              int L, int spin,
              ssht_dl_method_t dl_method,
              int verbosity)
        void ssht_core_mw_forward_sov_conv_sym_ss_real(double complex *flm, const double *f,
              int L,
              ssht_dl_method_t dl_method,
              int verbosity)
        void ssht_core_dh_inverse_sov(double complex *f, const double complex *flm,
                                  int L, int spin, int verbosity)
        void ssht_core_dh_inverse_sov_real(double *f, const double complex *flm,
                                  int L, int verbosity)
        void ssht_core_dh_forward_sov(double complex *flm, const double complex *f,
                                  int L, int spin, int verbosity)
        void ssht_core_dh_forward_sov_real(double complex *flm, const double *f,
                                  int L, int verbosity)
        void ssht_core_gl_inverse_sov(double complex *f, const double complex *flm,
                                  int L, int spin, int verbosity)
        void ssht_core_gl_inverse_sov_real(double *f, const double complex *flm,
                                  int L, int verbosity)
        void ssht_core_gl_forward_sov(double complex *flm, const double complex *f,
                                  int L, int spin, int verbosity)
        void ssht_core_gl_forward_sov_real(double complex *flm, const double *f,
                                  int L, int verbosity)
        double ssht_sampling_mw_t2theta(int t, int L)
        double ssht_sampling_mw_p2phi(int p, int L)
        double ssht_sampling_mw_ss_t2theta(int t, int L)
        double ssht_sampling_mw_ss_p2phi(int p, int L)
        double ssht_sampling_dh_t2theta(int t, int L)
        double ssht_sampling_dh_p2phi(int p, int L)
        void ssht_sampling_gl_thetas_weights(double *thetas, double *weights, int L)
        double ssht_sampling_gl_p2phi(int p, int L)

# adjoints

        void ssht_adjoint_mw_inverse_sov_sym(double complex *flm, 
             double complex *f, 
             int L, int spin, 
             ssht_dl_method_t dl_method,
             int verbosity)
        void ssht_adjoint_mw_inverse_sov_sym_real(double complex *flm, 
            double *f, 
            int L,
            ssht_dl_method_t dl_method, 
            int verbosity)
        void ssht_adjoint_mw_forward_sov_sym(double complex *f, 
             double complex *flm,
             int L, int spin,
             ssht_dl_method_t dl_method,
             int verbosity)
        void ssht_adjoint_mw_forward_sov_sym_real(double *f, 
            double complex *flm,
            int L,
            ssht_dl_method_t dl_method,
            int verbosity)

        void ssht_adjoint_mw_inverse_sov_sym_pole(double complex *flm, double complex *f,
            double complex f_sp, double phi_sp,
            int L, int spin, 
            ssht_dl_method_t dl_method,
            int verbosity)
        void ssht_adjoint_mw_inverse_sov_sym_real_pole(double complex *flm, 
                 double *f, 
                 double f_sp,
                 int L, 
                 ssht_dl_method_t dl_method,
                 int verbosity)
        void ssht_adjoint_mw_forward_sov_sym_pole(double complex *f, 
            double complex *f_sp, double *phi_sp,
            double complex *flm, 
            int L, int spin, 
            ssht_dl_method_t dl_method,
            int verbosity)
        void ssht_adjoint_mw_forward_sov_sym_real_pole(double *f, 
                 double *f_sp,
                 double complex *flm, 
                 int L, 
                 ssht_dl_method_t dl_method, 
                 int verbosity)
        void ssht_adjoint_mw_inverse_sov_sym_ss(double complex *flm, double complex *f, 
          int L, int spin, 
          ssht_dl_method_t dl_method,
          int verbosity)
        void ssht_adjoint_mw_inverse_sov_sym_ss_real(double complex *flm, double *f, 
               int L, 
               ssht_dl_method_t dl_method, 
               int verbosity)
        void ssht_adjoint_mw_forward_sov_sym_ss(double complex *f, double complex *flm,
          int L, int spin,
          ssht_dl_method_t dl_method,
          int verbosity)
        void ssht_adjoint_mw_forward_sov_sym_ss_real(double *f, 
               double complex *flm,
               int L,
               ssht_dl_method_t dl_method,
               int verbosity)



#----------------------------------------------------------------------------------------------------#

def ssht_forward_mw_complex(np.ndarray[ double complex, ndim=2, mode="c"] f_mw_c not None,int L,int spin):

        cdef ssht_dl_method_t dl_method = SSHT_DL_RISBO
        f_lm = np.empty([L * L,], dtype=complex)
        ssht_core_mw_forward_sov_conv_sym(<double complex*> np.PyArray_DATA(f_lm),<const double complex*> np.PyArray_DATA(f_mw_c), L, spin, dl_method, 0)
        return f_lm
        
def ssht_inverse_mw_complex(np.ndarray[ double complex, ndim=1, mode="c"] f_lm not None, int L, int spin):

        cdef ssht_dl_method_t dl_method = SSHT_DL_RISBO
        f_mw_c = np.empty([L,2*L-1,], dtype=complex)
        ssht_core_mw_inverse_sov_sym(<double complex*> np.PyArray_DATA(f_mw_c),<const double complex*> np.PyArray_DATA(f_lm), L, spin, dl_method, 0)
        return f_mw_c


def ssht_inverse_mw_complex_adjoint(np.ndarray[ double complex, ndim=2, mode="c"] f_mw_c not None,int L,int spin):

        cdef ssht_dl_method_t dl_method = SSHT_DL_RISBO
        f_lm = np.empty([L * L,], dtype=complex)
        ssht_adjoint_mw_inverse_sov_sym(<double complex*> np.PyArray_DATA(f_lm),<const double complex*> np.PyArray_DATA(f_mw_c), L, spin, dl_method, 0)
        return f_lm
        
def ssht_forward_mw_complex_adjoint(np.ndarray[ double complex, ndim=1, mode="c"] f_lm not None, int L, int spin):

        cdef ssht_dl_method_t dl_method = SSHT_DL_RISBO
        f_mw_c = np.empty([L,2*L-1,], dtype=complex)
        ssht_adjoint_mw_forward_sov_sym(<double complex*> np.PyArray_DATA(f_mw_c),<const double complex*> np.PyArray_DATA(f_lm), L, spin, dl_method, 0)
        return f_mw_c

#----------------------------------------------------------------------------------------------------#

def ssht_forward_mw_real(np.ndarray[ double, ndim=2, mode="c"] f_mw_r not None,int L):

        cdef ssht_dl_method_t dl_method = SSHT_DL_RISBO
        f_lm = np.empty([L * L,], dtype=complex)
        ssht_core_mw_forward_sov_conv_sym_real(<double complex*> np.PyArray_DATA(f_lm),<const double*> np.PyArray_DATA(f_mw_r), L, dl_method, 0)
        return f_lm
        
def ssht_inverse_mw_real(np.ndarray[ double complex, ndim=1, mode="c"] f_lm not None, int L):

        cdef ssht_dl_method_t dl_method = SSHT_DL_RISBO
        f_mw_r = np.empty([L,2*L-1,], dtype=np.float_)
        ssht_core_mw_inverse_sov_sym_real(<double*> np.PyArray_DATA(f_mw_r),<const double complex*> np.PyArray_DATA(f_lm), L, dl_method, 0)
        return f_mw_r

def ssht_inverse_mw_real_adjoint(np.ndarray[ double, ndim=2, mode="c"] f_mw_r not None,int L):

        cdef ssht_dl_method_t dl_method = SSHT_DL_RISBO
        f_lm = np.empty([L * L,], dtype=complex)
        ssht_adjoint_mw_inverse_sov_sym_real(<double complex*> np.PyArray_DATA(f_lm),<const double*> np.PyArray_DATA(f_mw_r), L, dl_method, 0)
        return f_lm
        
def ssht_forward_mw_real_adjoint(np.ndarray[ double complex, ndim=1, mode="c"] f_lm not None, int L):

        cdef ssht_dl_method_t dl_method = SSHT_DL_RISBO
        f_mw_r = np.empty([L,2*L-1,], dtype=np.float_)
        ssht_adjoint_mw_forward_sov_sym_real(<double*> np.PyArray_DATA(f_mw_r),<const double complex*> np.PyArray_DATA(f_lm), L, dl_method, 0)
        return f_mw_r

#----------------------------------------------------------------------------------------------------#

def ssht_forward_mwss_complex(np.ndarray[ double complex, ndim=2, mode="c"] f_mwss_c not None,int L,int spin):

        cdef ssht_dl_method_t dl_method = SSHT_DL_RISBO
        f_lm = np.empty([L * L,], dtype=complex)
        ssht_core_mw_forward_sov_conv_sym_ss(<double complex*> np.PyArray_DATA(f_lm),<const double complex*> np.PyArray_DATA(f_mwss_c), L, spin, dl_method, 0)
        return f_lm
        
def ssht_inverse_mwss_complex(np.ndarray[ double complex, ndim=1, mode="c"] f_lm not None, int L, int spin):

        cdef ssht_dl_method_t dl_method = SSHT_DL_RISBO
        f_mwss_c = np.empty([L+1,2*L,], dtype=complex)
        ssht_core_mw_inverse_sov_sym_ss(<double complex*> np.PyArray_DATA(f_mwss_c),<const double complex*> np.PyArray_DATA(f_lm), L, spin, dl_method, 0)
        return f_mwss_c


def ssht_inverse_mwss_complex_adjoint(np.ndarray[ double complex, ndim=2, mode="c"] f_mwss_c not None,int L,int spin):

        cdef ssht_dl_method_t dl_method = SSHT_DL_RISBO
        f_lm = np.empty([L * L,], dtype=complex)
        ssht_adjoint_mw_inverse_sov_sym_ss(<double complex*> np.PyArray_DATA(f_lm),<const double complex*> np.PyArray_DATA(f_mwss_c), L, spin, dl_method, 0)
        return f_lm
        
def ssht_forward_mwss_complex_adjoint(np.ndarray[ double complex, ndim=1, mode="c"] f_lm not None, int L, int spin):

        cdef ssht_dl_method_t dl_method = SSHT_DL_RISBO
        f_mwss_c = np.empty([L+1,2*L,], dtype=complex)
        ssht_adjoint_mw_forward_sov_sym_ss(<double complex*> np.PyArray_DATA(f_mwss_c),<const double complex*> np.PyArray_DATA(f_lm), L, spin, dl_method, 0)
        return f_mwss_c
#----------------------------------------------------------------------------------------------------#

def ssht_forward_mwss_real(np.ndarray[ double, ndim=2, mode="c"] f_mwss_r not None,int L):

        cdef ssht_dl_method_t dl_method = SSHT_DL_RISBO
        f_lm = np.empty([L * L,], dtype=complex)
        ssht_core_mw_forward_sov_conv_sym_ss_real(<double complex*> np.PyArray_DATA(f_lm),<const double*> np.PyArray_DATA(f_mwss_r), L, dl_method, 0)
        return f_lm
        
def ssht_inverse_mwss_real(np.ndarray[ double complex, ndim=1, mode="c"] f_lm not None, int L):

        cdef ssht_dl_method_t dl_method = SSHT_DL_RISBO
        f_mwss_r = np.empty([L+1,2*L,], dtype=np.float_)
        ssht_core_mw_inverse_sov_sym_ss_real(<double*> np.PyArray_DATA(f_mwss_r),<const double complex*> np.PyArray_DATA(f_lm), L, dl_method, 0)
        return f_mwss_r

def ssht_inverse_mwss_real_adjoint(np.ndarray[ double, ndim=2, mode="c"] f_mwss_r not None,int L):

        cdef ssht_dl_method_t dl_method = SSHT_DL_RISBO
        f_lm = np.empty([L * L,], dtype=complex)
        ssht_adjoint_mw_inverse_sov_sym_ss_real(<double complex*> np.PyArray_DATA(f_lm),<const double*> np.PyArray_DATA(f_mwss_r), L, dl_method, 0)
        return f_lm
        
def ssht_forward_mwss_real_adjoint(np.ndarray[ double complex, ndim=1, mode="c"] f_lm not None, int L):

        cdef ssht_dl_method_t dl_method = SSHT_DL_RISBO
        f_mwss_r = np.empty([L+1,2*L,], dtype=np.float_)
        ssht_adjoint_mw_forward_sov_sym_ss_real(<double*> np.PyArray_DATA(f_mwss_r),<const double complex*> np.PyArray_DATA(f_lm), L, dl_method, 0)
        return f_mwss_r

#----------------------------------------------------------------------------------------------------#

def ssht_forward_mw_complex_pole(np.ndarray[ double complex, ndim=2, mode="c"] f_mw_c not None, double complex f_sp, double phi_sp, int L,int spin):

        cdef ssht_dl_method_t dl_method = SSHT_DL_RISBO
        f_lm = np.empty([L * L,], dtype=complex)
        ssht_core_mw_forward_sov_conv_sym_pole(<double complex*> np.PyArray_DATA(f_lm),<const double complex*> np.PyArray_DATA(f_mw_c),  f_sp,  phi_sp, L, spin, dl_method, 0)
        return f_lm
        
def ssht_inverse_mw_complex_pole(np.ndarray[ double complex, ndim=1, mode="c"] f_lm not None, int L, int spin):

        cdef ssht_dl_method_t dl_method = SSHT_DL_RISBO
        f_mw_c = np.empty([L-1,(2*L-1),], dtype=complex)
        cdef double complex f_sp
        cdef double phi_sp
        ssht_core_mw_inverse_sov_sym_pole(<double complex*> np.PyArray_DATA(f_mw_c),  &f_sp,  &phi_sp, <const double complex*> np.PyArray_DATA(f_lm), L, spin, dl_method, 0)
        return f_mw_c, f_sp, phi_sp


def ssht_inverse_mw_complex_pole_adjoint(np.ndarray[ double complex, ndim=2, mode="c"] f_mw_c not None, double complex f_sp, double phi_sp, int L,int spin):

        cdef ssht_dl_method_t dl_method = SSHT_DL_RISBO
        f_lm = np.empty([L * L,], dtype=complex)
        ssht_adjoint_mw_inverse_sov_sym_pole(<double complex*> np.PyArray_DATA(f_lm),<const double complex*> np.PyArray_DATA(f_mw_c),  f_sp,  phi_sp, L, spin, dl_method, 0)
        return f_lm
        
def ssht_forward_mw_complex_pole_adjoint(np.ndarray[ double complex, ndim=1, mode="c"] f_lm not None, int L, int spin):

        cdef ssht_dl_method_t dl_method = SSHT_DL_RISBO
        f_mw_c = np.empty([L-1,(2*L-1),], dtype=complex)
        cdef double complex f_sp
        cdef double phi_sp
        ssht_adjoint_mw_forward_sov_sym_pole(<double complex*> np.PyArray_DATA(f_mw_c),  &f_sp,  &phi_sp, <const double complex*> np.PyArray_DATA(f_lm), L, spin, dl_method, 0)
        return f_mw_c, f_sp, phi_sp

#----------------------------------------------------------------------------------------------------#

def ssht_forward_mw_real_pole(np.ndarray[ double, ndim=2, mode="c"] f_mw_r not None, double f_sp,  int L):

        cdef ssht_dl_method_t dl_method = SSHT_DL_RISBO
        f_lm = np.empty([L * L,], dtype=complex)
        ssht_core_mw_forward_sov_conv_sym_real_pole(<double complex*> np.PyArray_DATA(f_lm),<const double*> np.PyArray_DATA(f_mw_r), f_sp, L, dl_method, 0)
        return f_lm
        
def ssht_inverse_mw_real_pole(np.ndarray[ double complex, ndim=1, mode="c"] f_lm not None, int L):

        cdef ssht_dl_method_t dl_method = SSHT_DL_RISBO
        f_mw_r = np.empty([L-1,2*L-1,], dtype=np.float_)
        cdef double f_sp
        ssht_core_mw_inverse_sov_sym_real_pole(<double*> np.PyArray_DATA(f_mw_r), &f_sp, <const double complex*> np.PyArray_DATA(f_lm), L, dl_method, 0)
        return f_mw_r, f_sp

def ssht_inverse_mw_real_pole_adjoint(np.ndarray[ double, ndim=2, mode="c"] f_mw_r not None, double f_sp,  int L):

        cdef ssht_dl_method_t dl_method = SSHT_DL_RISBO
        f_lm = np.empty([L * L,], dtype=complex)
        ssht_adjoint_mw_inverse_sov_sym_real_pole(<double complex*> np.PyArray_DATA(f_lm),<const double*> np.PyArray_DATA(f_mw_r), f_sp, L, dl_method, 0)
        return f_lm
        
def ssht_forward_mw_real_pole_adjoint(np.ndarray[ double complex, ndim=1, mode="c"] f_lm not None, int L):

        cdef ssht_dl_method_t dl_method = SSHT_DL_RISBO
        f_mw_r = np.empty([L-1,2*L-1,], dtype=np.float_)
        cdef double f_sp
        ssht_adjoint_mw_forward_sov_sym_real_pole(<double*> np.PyArray_DATA(f_mw_r), &f_sp, <const double complex*> np.PyArray_DATA(f_lm), L, dl_method, 0)
        return f_mw_r, f_sp

#----------------------------------------------------------------------------------------------------#

def ssht_forward_dh_complex(np.ndarray[ double complex, ndim=2, mode="c"] f_dh_c not None,int L,int spin):

        f_lm = np.empty([L * L,], dtype=complex)
        ssht_core_dh_forward_sov(<double complex*> np.PyArray_DATA(f_lm),<const double complex*> np.PyArray_DATA(f_dh_c), L, spin, 0)
        return f_lm
        
def ssht_inverse_dh_complex(np.ndarray[ double complex, ndim=1, mode="c"] f_lm not None, int L, int spin):

        f_dh_c = np.empty([2*L,2*L-1,], dtype=complex)
        ssht_core_dh_inverse_sov(<double complex*> np.PyArray_DATA(f_dh_c),<const double complex*> np.PyArray_DATA(f_lm), L, spin, 0)
        return f_dh_c

#----------------------------------------------------------------------------------------------------#

def ssht_forward_dh_real(np.ndarray[ double, ndim=2, mode="c"] f_dh_r not None,int L):

        f_lm = np.empty([L * L,], dtype=complex)
        ssht_core_dh_forward_sov_real(<double complex*> np.PyArray_DATA(f_lm),<const double*> np.PyArray_DATA(f_dh_r), L,  0)
        return f_lm
        
def ssht_inverse_dh_real(np.ndarray[ double complex, ndim=1, mode="c"] f_lm not None, int L):

        f_dh_r = np.empty([2*L,2*L-1,], dtype=np.float_)
        ssht_core_dh_inverse_sov_real(<double*> np.PyArray_DATA(f_dh_r),<const double complex*> np.PyArray_DATA(f_lm), L, 0)
        return f_dh_r


#----------------------------------------------------------------------------------------------------#

def ssht_forward_gl_complex(np.ndarray[ double complex, ndim=2, mode="c"] f_gl_c not None,int L,int spin):

        cdef ssht_dl_method_t dl_method = SSHT_DL_RISBO
        f_lm = np.empty([L * L,], dtype=complex)
        ssht_core_gl_forward_sov(<double complex*> np.PyArray_DATA(f_lm),<const double complex*> np.PyArray_DATA(f_gl_c), L, spin, 0)
        return f_lm
        
def ssht_inverse_gl_complex(np.ndarray[ double complex, ndim=1, mode="c"] f_lm not None, int L, int spin):

        cdef ssht_dl_method_t dl_method = SSHT_DL_RISBO
        f_gl_c = np.empty([L,2*L-1,], dtype=complex)
        ssht_core_gl_inverse_sov(<double complex*> np.PyArray_DATA(f_gl_c),<const double complex*> np.PyArray_DATA(f_lm), L, spin, 0)
        return f_gl_c

#----------------------------------------------------------------------------------------------------#

def ssht_forward_gl_real(np.ndarray[ double, ndim=2, mode="c"] f_gl_r not None,int L):

        cdef ssht_dl_method_t dl_method = SSHT_DL_RISBO
        f_lm = np.empty([L * L,], dtype=complex)
        ssht_core_gl_forward_sov_real(<double complex*> np.PyArray_DATA(f_lm),<const double*> np.PyArray_DATA(f_gl_r), L, 0)
        return f_lm
        
def ssht_inverse_gl_real(np.ndarray[ double complex, ndim=1, mode="c"] f_lm not None, int L):

        cdef ssht_dl_method_t dl_method = SSHT_DL_RISBO
        f_gl_r = np.empty([L,2*L-1,], dtype=np.float_)
        ssht_core_gl_inverse_sov_real(<double*> np.PyArray_DATA(f_gl_r),<const double complex*> np.PyArray_DATA(f_lm), L, 0)
        return f_gl_r



#----------------------------------------------------------------------------------------------------#

# easy to use functions to perform forward and backward transfroms

def forward(f, int L, int Spin=0, str Method='MW', bint Reality=False, str backend="SSHT", **kwargs):
    from pyssht.ducc_interface import forward as _ducc0_forward
    from pyssht.parameters import method, Ducc

    params = method(Method, spin=Spin, reality=Reality, backend=backend, **kwargs)
    if params.method == 'MW_POLE':
        if Reality:
            f, f_sp = f
        else:
            f, f_sp, phi_sp = f

    if f.ndim != 2:
      raise ssht_input_error('f must be 2D numpy array')

    if isinstance(params, Ducc):
        return _ducc0_forward(f, L, **params.asDict())

    if f.dtype == np.float_ and params.reality == False:
        print('Real signal given but Reality flag is False. Set Reality = True to improve performance')
        f_new = np.empty(sample_shape(L,Method=Method), dtype=np.complex_)
        f_new = f + 1j*np.zeros(sample_shape(L,Method=Method), dtype=np.float_)
        f = f_new

    if f.dtype == np.complex_ and params.reality == True:
        print('Complex signal given but Reality flag is True. Ignoring complex component')
        f_new = np.real(f)
        f = f_new.copy(order='c')

    # do correct transform
    if params.method == 'MW':
        if params.reality:
            flm = ssht_forward_mw_real(f,L)
        else:
            flm = ssht_forward_mw_complex(f,L,params.spin)
            
    if params.method == 'MW_POLE':
        if params.reality:
            flm = ssht_forward_mw_real_pole(f,f_sp,L)
        else:
            flm = ssht_forward_mw_complex_pole(f,f_sp,phi_sp,L,params.spin)
            
    if params.method == 'MWSS':
        if params.reality:
            flm = ssht_forward_mwss_real(f,L)
        else:
            flm = ssht_forward_mwss_complex(f,L,params.spin)

    if params.method == 'DH':
        if params.reality:
            flm = ssht_forward_dh_real(f,L)
        else:
            flm = ssht_forward_dh_complex(f,L,params.spin)
            
    if params.method == 'GL':
        if params.reality:
            flm = ssht_forward_gl_real(f,L)
        else:
            flm = ssht_forward_gl_complex(f,L,params.spin)
            
            
    return flm


def inverse(flm, int L, int Spin=0, str Method='MW', bint Reality=False, str backend="SSHT", **kwargs):
    from pyssht.ducc_interface import inverse as _ducc0_inverse
    from pyssht.parameters import method, Ducc

    if flm.ndim != 1:
      raise ssht_input_error('flm must be 1D numpy array')

    params = method(Method, spin=Spin, reality=Reality, backend=backend, **kwargs)
    
    if isinstance(params, Ducc):
        return _ducc0_inverse(flm, L, params.spin, params.method, params.reality, params.nthreads)

    # do correct transform
    if params.method == 'MW':
        if params.reality:
            f = ssht_inverse_mw_real(flm,L)
        else:
            f = ssht_inverse_mw_complex(flm,L,params.spin)
            

    if params.method == 'MW_POLE':
        if params.reality:
            f = ssht_inverse_mw_real_pole(flm,L)
        else:
            f = ssht_inverse_mw_complex_pole(flm,L,params.spin)

    if params.method == 'MWSS':
        if params.reality:
            f = ssht_inverse_mwss_real(flm,L)
        else:
            f = ssht_inverse_mwss_complex(flm,L,params.spin)
          
    if params.method == 'DH':
        if params.reality:
            f = ssht_inverse_dh_real(flm,L)
        else:
            f = ssht_inverse_dh_complex(flm,L,params.spin)
            
    if params.method == 'GL':
        if params.reality:
            f = ssht_inverse_gl_real(flm,L)
        else:
            f = ssht_inverse_gl_complex(flm,L,params.spin)
            
            
            
    return f

def inverse_adjoint(f, int L, int Spin=0, str Method='MW', bint Reality=False, str backend="SSHT", **kwargs):
    from pyssht.ducc_interface import inverse_adjoint as _ducc0_inverse_adjoint
    from pyssht.parameters import method, Ducc

    params = method(Method, spin=Spin, reality=Reality, backend=backend)

    if isinstance(params, Ducc):
        return _ducc0_inverse_adjoint(f, L, params.spin, params.method, params.reality, params.nthreads)

    if params.method == 'MW_POLE':
        if Reality:
            f, f_sp = f
        else:
            f, f_sp, phi_sp = f

    if f.ndim != 2:
      raise ssht_input_error('f must be 2D numpy array')

    if params.method not in ('MW', 'MWSS'):
        raise ssht_input_error('Method is not recognised, Methods are: MW, MWSS')

    if f.dtype == np.float_ and params.reality == False:
        print('Real signal given but Reality flag is False. Set Reality = True to improve performance')
        f_new = np.empty(sample_shape(L,Method=params.method), dtype=np.complex_)
        f_new = f + 1j*np.zeros(sample_shape(L,Method=params.method), dtype=np.float_)
        f = f_new

    if f.dtype == np.complex_ and params.reality == True:
        print('Complex signal given but Reality flag is True. Ignoring complex component')
        f_new = np.real(f)
        f = f_new.copy(order='c')
        
    # do correct transform
    if params.method == 'MW':
        if params.reality:
            flm = ssht_inverse_mw_real_adjoint(f,L)
        else:
            flm = ssht_inverse_mw_complex_adjoint(f,L,params.spin)
    elif params.method == 'MWSS':
        if params.reality:
            flm = ssht_inverse_mwss_real_adjoint(f,L)
        else:
            flm = ssht_inverse_mwss_complex_adjoint(f,L,params.spin)
            
            
    return flm


def forward_adjoint(flm, int L, int Spin=0, str Method='MW', bint Reality=False, str backend="SSHT", **kwargs):
    from pyssht.ducc_interface import forward_adjoint as _ducc0_forward_adjoint
    from pyssht.parameters import method, Ducc

    params = method(Method, spin=Spin, reality=Reality, backend=backend)

    if isinstance(params, Ducc):
        return _ducc0_forward_adjoint(flm, L, params.spin, params.method, params.reality, params.nthreads)

    if flm.ndim != 1:
      raise ssht_input_error('flm must be 1D numpy array')

    if params.method not in ( 'MW', 'MWSS'):
        raise ssht_input_error('Method is not recognised, Methods are: MW, MWSS')
    
    if params.method == 'MW':
        if params.reality:
            f = ssht_forward_mw_real_adjoint(flm,L)
        else:
            f = ssht_forward_mw_complex_adjoint(flm,L,params.spin)
    elif params.method == 'MWSS':
         if params.reality:
             f = ssht_forward_mwss_real_adjoint(flm,L)
         else:
             f = ssht_forward_mwss_complex_adjoint(flm,L,params.spin)
                      
    return f
#----------------------------------------------------------------------------------------------------#

# index to ell em and back function

def isqrt(int n):
    cdef int square = 1, delta = 3
    while square <= n:
        square += delta
        delta  += 2
    return (delta/2 -1)

def elm2ind( int el, int m):

  return el * el + el + m


def ind2elm(int ind):

  cdef int ell, em
  el = cy_isqrt(ind)
  em = ind - (el)*(el) - (el)

  return el, em

# index to ell em and back function

cdef inline int cy_isqrt(int n):
    cdef int square = 1, delta = 3
    while square <= n:
        square += delta
        delta  += 2
    return (int(delta/2) -1)

cdef inline int cy_elm2ind( int el, int m):

  return el * el + el + m


def theta_to_index(double theta, int L, str Method="MW"):
  cdef int p
  cdef np.ndarray[np.float_t, ndim=1] theta_gl_grid, phi_gl_grid

  if Method == 'GL':
    theta_gl_grid, phi_gl_grid = sample_positions(L,Method="GL")

  if Method == 'MW':
    p = int((theta*(2*L-1)/pi-1)/2)  # (2.0*t + 1.0) * SSHT_PI / (2.0*L - 1.0)
  if Method == 'MWSS':
    p = int((theta*(2*L)/pi)/2)  # 2.0 * t * SSHT_PI / (2.0 * L)
  if Method == 'DH':
    p = int((theta*(4*L)/pi-1)/2)  # (2.0*t + 1.0) * SSHT_PI / (4.0*L)
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
    q = int(phi*(2*L-1)/(2*pi))      # 2.0 * p * SSHT_PI / (2.0*L - 1.0)
  if Method == 'MWSS':
    q = int(phi*(2*L)/(2*pi))      # 2.0 * p * SSHT_PI / (2.0*L)
  if Method == 'DH':
    q = int(phi*(2*L-1)/(2*pi))      # 2.0 * p * SSHT_PI / (2.0*L - 1.0)
  if Method == 'GL':
    q = int(phi*(2*L-1)/(2*pi))      # 2.0 * p * SSHT_PI / (2.0*L - 1.0)
  return q

cdef inline int cy_theta_to_index(double theta, int L, METHOD_TYPE Method_enum):
  cdef int p

  if Method_enum == GL:
    theta_gl_grid, phi_gl_grid = sample_positions(L,Method="GL")

  if Method_enum == MW:
    p = int((theta*(2*L-1)/pi-1)/2+0.5)  # (2.0*t + 1.0) * SSHT_PI / (2.0*L - 1.0)
  if Method_enum == MWSS:
    p = int((theta*(2*L)/pi)/2+0.5)  # 2.0 * t * SSHT_PI / (2.0 * L)
  if Method_enum == DH:
    p = int((theta*(4*L)/pi-1)/2+0.5)  # (2.0*t + 1.0) * SSHT_PI / (4.0*L)
  if Method_enum == GL:
    if theta > theta_gl_grid[L-1]:
      p = L-1
    for k in range(theta_gl_grid.size):
      if theta < theta_gl_grid[k]:
        p = k-1
        break
  return p

def index_to_theta(int index, int L, str Method):
  cdef METHOD_TYPE Method_enum=get_method_enum(Method)
  return cy_index_to_theta(index, L, Method_enum)

cdef inline double cy_index_to_theta(int index, int L, METHOD_TYPE Method_enum):
  cdef double p

  if Method_enum == MW:
    p = pi * float(2*index+1) / ( 2.0 * float(L) - 1.0 )# (2.0*t + 1.0) * SSHT_PI / (2.0*L - 1.0)
  if Method_enum == MWSS: 
    p = pi*float(index) / float(L) # 2.0 * t * SSHT_PI / (2.0 * L);

  return p

cdef inline int cy_phi_to_index(double phi, int L, METHOD_TYPE Method_enum):
  cdef int q
  
  if Method_enum == MW:
    q = int(phi*(2*L-1)/(2*pi)+0.5)      # 2.0 * p * SSHT_PI / (2.0*L - 1.0)
    if q == 2*L-1:
      q = 0
  if Method_enum == MWSS:
    q = int(phi*(2*L)/(2*pi)+0.5)      # 2.0 * p * SSHT_PI / (2.0*L)
    if q == 2*L:
      q = 0
  if Method_enum == DH:
    q = int(phi*(2*L-1)/(2*pi)+0.5)      # 2.0 * p * SSHT_PI / (2.0*L - 1.0)
    if q == 2*L-1:
      q = 0
  if Method_enum == GL:
    q = int(phi*(2*L-1)/(2*pi)+0.5)      # 2.0 * p * SSHT_PI / (2.0*L - 1.0)
    if q == 2*L-1:
      q = 0

  return q

def index_to_phi(int index, int L, str Method):
  cdef METHOD_TYPE Method_enum=get_method_enum(Method)
  return cy_index_to_phi(index, L, Method_enum)


cdef inline double cy_index_to_phi(int index, int L, METHOD_TYPE Method_enum):
  cdef double q

  if Method_enum == MW:
    q = 2.0 * pi * float(index) / ( 2.0 * float(L) - 1.0 )
  if Method_enum == MWSS:
    q = pi * float(index) / float(L)

  return q


def sample_length(int L, Method = 'MW'):
  from pyssht.parameters import method
  params = method(Method)
  cdef int n_theta, n_phi
  n_theta, n_phi = sample_shape(L, Method=params.method)

  return n_theta*n_phi
# get the shape of the signal on the sphere for different sampling theorems

def sample_shape(int L, str Method = 'MW'):
  from pyssht.parameters import method
  params = method(Method)

  cdef int n_theta, n_phi
  if params.method == 'MW':
    n_theta = L
    n_phi   = 2*L-1

  if params.method == 'MW_POLE':
    n_theta = L-1
    n_phi   = 2*L-1

  if params.method == 'MWSS':
    n_theta = L+1
    n_phi   = 2*L

  if params.method == 'GL':
    n_theta = L
    n_phi   = 2*L-1

  if params.method == 'DH':
    n_theta = 2*L
    n_phi   = 2*L-1



  return (n_theta, n_phi)

def sample_positions(int L, str Method = 'MW', bint Grid=False):
  from pyssht.parameters import method
  params = method(Method)

  cdef int n_theta, n_phi, i
  cdef np.ndarray[np.float_t, ndim=1] thetas, phis

  n_theta, n_phi = sample_shape(L, Method=params.method)
  thetas = np.empty(n_theta, dtype=np.float_)
  phis   = np.empty(n_phi,   dtype=np.float_)


  if params.method == 'MW':
    for i in range(n_theta):
      thetas[i] = ssht_sampling_mw_t2theta(i, L)
    for i in range(n_phi):
      phis[i] = ssht_sampling_mw_p2phi(i, L)

  if params.method == 'MWSS':
    for i in range(n_theta):
      thetas[i] = ssht_sampling_mw_ss_t2theta(i, L)
    for i in range(n_phi):
      phis[i] = ssht_sampling_mw_ss_p2phi(i, L)

  if params.method == 'DH':
    for i in range(n_theta):
      thetas[i] = ssht_sampling_dh_t2theta(i, L)
    for i in range(n_phi):
      phis[i] = ssht_sampling_dh_p2phi(i, L)

  if params.method == 'GL':
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
    phi[phi<0] += 2*pi 

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
    phi[phi<0] += 2*pi 

  return (r, theta, phi)

def theta_phi_to_ra_dec(theta, phi, bint Degrees=False):
  dec = (theta - pi/2)*(-1)
  ra  = phi

  if Degrees:
    dec = dec*180/pi
    ra  = ra*180/pi

  return dec, ra


def ra_dec_to_theta_phi(ra, dec, bint Degrees=False):
  if Degrees:
    dec = dec*pi/180
    ra  = ra*pi/180

  theta = pi/2 - dec
  phi  = ra


  return theta, phi

def old_make_rotation_matrix(list rot):
  cdef np.ndarray[np.float_t, ndim=2] rot_matix = np.empty((3,3), dtype=np.float_)

  cdef double c1 = cos(float(rot[0]))
  cdef double s1 = sin(float(rot[0]))
  cdef double c2 = cos(float(rot[1]))
  cdef double s2 = sin(float(rot[1]))
  cdef double c3 = cos(float(rot[2]))
  cdef double s3 = sin(float(rot[2]))

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


def make_rotation_matrix(list rot):
  '''
  ZYZ rotation: R = Rz(alpha)Ry(beta)Rz(gamma)
  list rot = [alpha,beta,gamma]
  '''
  cdef np.ndarray[np.float_t, ndim=2] rot_matix = np.empty((3,3), dtype=np.float_)

  cdef double c1 = cos(float(rot[0]))
  cdef double s1 = sin(float(rot[0]))
  cdef double c2 = cos(float(rot[1]))
  cdef double s2 = sin(float(rot[1]))
  cdef double c3 = cos(float(rot[2]))
  cdef double s3 = sin(float(rot[2]))

  rot_matix[0,0] = c1*c2*c3 - s1*s3
  rot_matix[0,1] = -c1*c2*s3 - c3*s1
  rot_matix[0,2] = c1*s2
  rot_matix[1,0] = c2*c3*s1 + c1*s3
  rot_matix[1,1] = c1*c3 - c2*s1*s3
  rot_matix[1,2] = s1*s2
  rot_matix[2,0] = -c3*s2
  rot_matix[2,1] = s2*s3
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


def rotate_image(image, rot_list, Method=None):
  if Method is not None:
    raise NotImplementedError()  

  Height = image.shape[0]
  Width = image.shape[1]
  # get theta, phi of pix positions
  i_array = np.arange(Height)
  j_array = np.arange(Width)

  j_array, i_array = np.meshgrid(j_array, i_array)

  theta = (i_array)*np.pi/Height
  phi = (j_array)*2*np.pi/Width

  # convert theta phi to x, y, z
  x, y, z = s2_to_cart(theta, phi)

  # rotate x, y, z
  x, y, z = rot_cart_2d(x, y, z, [-u for u in rot_list[::-1]])

  # convert x, y, z, to theta, phi
  theta, phi = cart_to_s2(x, y, z)

  # convert to indexes
  i_array = theta*(Height)/np.pi 
  j_array = phi*(Width)/(2*np.pi)

  # interpolate array
  new_image = bilinear_interpolate(image, j_array, i_array)


  return new_image

def bilinear_interpolate(im, x, y):
  x = np.asarray(x)
  y = np.asarray(y)

  x0 = np.floor(x).astype(int)
  x1 = x0 + 1
  y0 = np.floor(y).astype(int)
  y1 = y0 + 1

  wa = (x1-x) * (y1-y)
  wb = (x1-x) * (y-y0)
  wc = (x-x0) * (y1-y)
  wd = (x-x0) * (y-y0)

  x0 = np.clip(x0, 0, im.shape[1]-1)
  x1 = np.clip(x1, 0, im.shape[1]-1)
  y0 = np.clip(y0, 0, im.shape[0]-1)
  y1 = np.clip(y1, 0, im.shape[0]-1)

  return wa*im[y0, x0] + wb*im[y1, x0] + wc*im[y0, x1] + wd*im[y1, x1]


#----------------------------------------------------------------------------------------------------#

# Plotting functions

def plot_sphere(f, int L, str Method='MW', bint Close=True, bint Parametric=False, list Parametric_Scaling=[0.0,0.5], \
                     Output_File=None, bint Show=True, bint Color_Bar=True, Units=None, Color_Range=None, \
                     int Axis=True): 
    import matplotlib.pyplot as plt
    from matplotlib import cm, colors, colorbar, gridspec
    from mpl_toolkits.mplot3d import Axes3D
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    from pyssht.parameters import method

    # add ability to choose color bar min max
    # and sort out shapes of the plots
    
    params = method(Method)
    if params.method == 'MW_POLE':
        if len(f) == 2:
            f, f_sp = f
        else:
            f, f_sp, phi_sp = f

    (thetas, phis) = sample_positions(L, Method=params.method, Grid=True)

    if (thetas.size != f.size):
        ssht_input_error('Band limit L deos not match that of f')

    f_plot = f.copy()

    f_max = f_plot.max()
    f_min = f_plot.min()

    if Color_Range is None:
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
        f_normalised = (f_plot - vmin/(vmax - vmin))*Parametric_Scaling[1]+Parametric_Scaling[0]


    # % Close plot.
    if Close:
        (n_theta,n_phi) = sample_shape(L,Method=params.method)
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
  dec = np.arcsin((2*aux_angle + np.sin(2*aux_angle))/pi)
  ra = pi*x/(2*np.sqrt(2.0)*np.cos(aux_angle))


  theta, phi = ra_dec_to_theta_phi(ra, dec)

  for i from 0 <= i < Ny:
    for j from 0 <= j < Nx:
      if phi[i,j]<-pi or phi[i,j]>pi:
        phi[i,j] = np.nan
      else:
        if phi[i,j] < 0.0:
          phi[i,j] += 2*pi

 
  if not(rot is None):
    if len(rot) == 1:
      for i from 0 <= i < Ny:
        for j from 0 <= j < Nx:
          if not(np.isnan(phi[i,j])):
            phi[i,j] += rot_angle
            if phi[i,j] > 2*pi:
              phi[i,j] -= 2*pi
            if phi[i,j] < 0:
              phi[i,j] += 2*pi

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

cdef METHOD_TYPE get_method_enum(str Method):
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
  return Method_enum



def mollweide_projection(f, int L, int resolution=500, rot=None,
                        zoom_region=[np.sqrt(2.0)*2,np.sqrt(2.0)], str Method="MW"):

  cdef int i, j, n_phi, n_theta
  cdef np.ndarray[np.float_t, ndim=2] f_real, f_imag
  cdef METHOD_TYPE Method_enum=get_method_enum(Method)

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

def equatorial_projection(f, int L, int resolution=500, rot=None,\
                        list zoom_region=[-1,-1], str Method="MW", str Projection="MERCATOR", int Spin=0):

  if not len(zoom_region)==2:
    raise ssht_input_error('zoom_region must be a python list of length 2')

  cdef METHOD_TYPE Method_enum
  cdef EQUATORIAL_PROJECTION_TYPE Equatorial_Projection_enum

  if Projection == "MERCATOR":
    Equatorial_Projection_enum = MERCATOR
    if zoom_region[0] < -0.5:
      zoom_region = [pi-1E-5,7*pi/16]
    if zoom_region[1] > pi/2-1E-5:
      zoom_region[1] = pi/2-1E-5
  elif Projection == "SINE":
    Equatorial_Projection_enum = SINE
    if zoom_region[0] < -0.5:
      zoom_region = [pi-1E-5,pi/2]
  else:
    ssht_input_error('Projection type not recognised, Projectios are MERCATOR and SINE')
  if zoom_region[0] > pi-1E-5:
    zoom_region[0] = pi-1E-5

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
    return equatorial_projection_work(f, L, resolution=resolution, rot=rot,\
                        zoom_region=zoom_region, Method_enum=Method_enum, Equatorial_Projection_enum=Equatorial_Projection_enum)
  elif f.dtype == complex:
    n_theta, n_phi = sample_shape(L,Method=Method)
    f_real = np.empty((n_theta,n_phi), dtype=np.float_)
    f_imag = np.empty((n_theta,n_phi), dtype=np.float_)
    for i in range(n_theta):
      for j in range(n_phi):
        f_real[i,j] = f[i,j].real
        f_imag[i,j] = f[i,j].imag

    f_real_plot, mask_real = equatorial_projection_work(f_real, L, resolution=resolution, rot=rot,\
                        zoom_region=zoom_region, Method_enum=Method_enum, Equatorial_Projection_enum=Equatorial_Projection_enum)
    f_imag_plot, mask_imag = equatorial_projection_work(f_imag, L, resolution=resolution, rot=rot,\
                        zoom_region=zoom_region, Method_enum=Method_enum, Equatorial_Projection_enum=Equatorial_Projection_enum)
    if Spin != 0:
      f_real_plot_keep = f_real_plot

      rotation_angle = equatorial_projection_rotation_angle(resolution, zoom_region, rot=rot, \
        Equatorial_Projection_enum=Equatorial_Projection_enum)

      f_real_plot = f_real_plot_keep*np.cos(Spin*rotation_angle) - f_imag_plot*np.sin(Spin*rotation_angle)
      f_imag_plot = f_imag_plot*np.cos(Spin*rotation_angle)      + f_real_plot_keep*np.sin(Spin*rotation_angle)

    return f_real_plot, mask_real, f_imag_plot, mask_imag
  else:
    raise ssht_input_error("f dtype must be float or complex")

def equatorial_projection_rotation_angle(int resolution, list zoom_region, rot=None, \
  EQUATORIAL_PROJECTION_TYPE Equatorial_Projection_enum=MERCATOR):
  if rot is None or len(rot) > 1:
    dummy = 0.0
  else:
    dummy = <double> rot[0]


  cdef double rot_angle = dummy

  cdef np.ndarray[np.float_t, ndim=2] rotation_angle, rot_matix
 
  cdef float x_pos, y_pos, z_pos, x_p_pos, y_p_pos, z_p_pos, rho, theta_pos, phi_pos, half_box_len_x, half_box_len_y, max_len
  cdef float delta_x, delta_y, delta_z, delta_theta, delta_phi, delta_y_plane, delta_x_plane
  cdef int n_theta, n_phi, n_theta_north, n_theta_south, i, j, p, q, i_rot, j_rot
  cdef int Nx, Ny

  half_box_len_x = forward_equatorial_projection_function_float_x(pi/2, zoom_region[0], Equatorial_Projection_enum)
  half_box_len_y = forward_equatorial_projection_function_float_y(pi/2-zoom_region[1], 0.0, Equatorial_Projection_enum)
 
  Nx = resolution
  Ny = int(<float>resolution*half_box_len_y/half_box_len_x)

  rotation_angle = np.full((Ny,Nx), np.nan, dtype=float)

  if not(rot is None):
    if len(rot) == 1:
      rot_list = [-rot[0],0.0,0.0]
    elif len(rot) == 3:
      rot_list = rot
    rot_matix = make_rotation_matrix(rot_list)
    delta_x = rot_matix[0,2]
    delta_y = rot_matix[1,2]
    delta_z = rot_matix[2,2]
  else:
    if Equatorial_Projection_enum == MERCATOR:
      return np.zeros((Ny,Nx), dtype=float)
    delta_x = 0.0
    delta_y = 0.0
    delta_z = 1.0

  for i in range(Nx):
    for j in range(Ny):
        x_pos     = (2.*half_box_len_x*(<float>i+0.5)/<float>Nx -half_box_len_x)
        y_pos     = (2.*half_box_len_y*(<float>j+0.5)/<float>Ny -half_box_len_y)
        theta_pos = inverse_equatorial_projection_function_float_theta(x_pos, y_pos, Equatorial_Projection_enum)
        phi_pos   = inverse_equatorial_projection_function_float_phi(x_pos, y_pos, Equatorial_Projection_enum)
        # print(i, j, x_pos, y_pos, rho, theta_pos)

        if not np.isnan(phi_pos):
          delta_theta  = sin(phi_pos)*delta_y+cos(phi_pos)*delta_x - tan(theta_pos)*delta_z
          delta_theta *= cos(theta_pos)

          delta_phi  = -tan(phi_pos)*delta_x + delta_y
          delta_phi *= cos(phi_pos)/sin(theta_pos)

          delta_x_plane  = forward_equatorial_projection_function_float_x_dtheta(theta_pos, phi_pos, \
            Equatorial_Projection_enum=Equatorial_Projection_enum)*delta_theta
          delta_x_plane += forward_equatorial_projection_function_float_x_dphi(theta_pos, phi_pos, \
            Equatorial_Projection_enum=Equatorial_Projection_enum)*delta_phi

          delta_y_plane  = forward_equatorial_projection_function_float_y_dtheta(theta_pos, phi_pos, \
            Equatorial_Projection_enum=Equatorial_Projection_enum)*delta_theta
          delta_y_plane += forward_equatorial_projection_function_float_y_dphi(theta_pos, phi_pos, \
            Equatorial_Projection_enum=Equatorial_Projection_enum)*delta_phi


          rotation_angle[j,i] = -atan2(delta_x_plane, delta_y_plane)


  return rotation_angle

def equatorial_projection_angle_array(int resolution, list zoom_region=[-1.,-1.], rot=None, str Projection="MERCATOR"):

  if rot is None:
    dummy = 0.0
  else:
    if len(rot) == 1:
      dummy = <double> rot[0]
    else:
      raise ssht_input_error("Three angle rotations are not supported for this function yet")
  cdef EQUATORIAL_PROJECTION_TYPE Equatorial_Projection_enum

  if Projection == "MERCATOR":
    Equatorial_Projection_enum = MERCATOR
    if zoom_region[0] < -0.5:
      zoom_region = [pi-1E-5,7*pi/16]
    if zoom_region[1] > pi/2-1E-5:
      zoom_region[1] = pi/2-1E-5
  elif Projection == "SINE":
    Equatorial_Projection_enum = SINE
    if zoom_region[0] < -0.5:
      zoom_region = [pi-1E-5,pi/2]
  else:
    ssht_input_error('Projection type not recognised, Projectios are MERCATOR and SINE')
  if zoom_region[0] > pi-1E-5:
    zoom_region[0] = pi-1E-5

  cdef double rot_angle = dummy

  cdef np.ndarray[np.float_t, ndim=2] theta_angle, phi_angle, rot_matix
 
  cdef float x_pos, y_pos, z_pos, x_p_pos, y_p_pos, z_p_pos, rho, theta_pos, phi_pos, half_box_len_x, half_box_len_y, max_len
  cdef int n_theta, n_phi, n_theta_north, n_theta_south, i, j, p, q, i_rot, j_rot
  cdef int Nx, Ny

  half_box_len_x = forward_equatorial_projection_function_float_x(pi/2, zoom_region[0], Equatorial_Projection_enum)
  half_box_len_y = forward_equatorial_projection_function_float_y(pi/2-zoom_region[1], 0.0, Equatorial_Projection_enum)
 
  Nx = resolution
  Ny = int(<float>resolution*half_box_len_y/half_box_len_x)

  
  theta_angle = np.full((Ny,Nx), np.nan, dtype=float)
  phi_angle   = np.full((Ny,Nx), np.nan, dtype=float)

  for i in range(Nx):
    for j in range(Ny):
        x_pos     = (2.*half_box_len_x*(<float>i+0.5)/<float>Nx -half_box_len_x)
        y_pos     = (2.*half_box_len_y*(<float>j+0.5)/<float>Ny -half_box_len_y)
        theta_pos = inverse_equatorial_projection_function_float_theta(x_pos, y_pos, Equatorial_Projection_enum)
        phi_pos   = inverse_equatorial_projection_function_float_phi(x_pos, y_pos, Equatorial_Projection_enum)
        # print(i, j, x_pos, y_pos, rho, theta_pos)

        if not np.isnan(phi_pos):
          phi_pos += rot_angle
          if phi_pos > 2*pi:
            phi_pos -= 2*pi
          if phi_pos < 0.0:
            phi_pos += 2*pi          
          
          theta_angle[j,i] = theta_pos
          phi_angle[j,i]   = phi_pos


  return theta_angle, phi_angle

def equatorial_projection_work(np.ndarray[ double, ndim=2, mode="c"] f, int L, int resolution=500, rot=None,\
                        list zoom_region=[-1,-1], METHOD_TYPE Method_enum=MW, \
                        EQUATORIAL_PROJECTION_TYPE Equatorial_Projection_enum=MERCATOR):


  if rot is None or len(rot) > 1:
    dummy = 0.0
  else:
    dummy = <double> rot[0]


  cdef double rot_angle = dummy

  cdef np.ndarray[np.float_t, ndim=2] theta, phi
  cdef np.ndarray[np.float_t, ndim=2] x_project, y_project
  cdef np.ndarray[np.float_t, ndim=2] xx, yy, zz, xx_p, yy_p, zz_p
  cdef np.ndarray[np.float_t, ndim=2] f_project, mask
  cdef np.ndarray[np.float_t, ndim=2] rot_matix

  cdef list rot_list

  cdef np.ndarray[np.int_t,   ndim=2] n_points

  cdef float x_pos, y_pos, z_pos, x_p_pos, y_p_pos, z_p_pos, rho, theta_pos, phi_pos, half_box_len_x, half_box_len_y, max_len
  
  cdef int n_theta, n_phi, n_theta_north, n_theta_south, i, j, p, q, i_rot, j_rot

  cdef int Nx, Ny

  cdef str Method

  if Method_enum==MW:
    Method="MW"
  elif Method_enum==MWSS:
    Method = "MWSS"
  elif Method_enum==DH:
    Method = "DH"
  elif Method_enum==GL:
    Method = "GL"

  half_box_len_x = forward_equatorial_projection_function_float_x(pi/2, zoom_region[0], Equatorial_Projection_enum)
  half_box_len_y = forward_equatorial_projection_function_float_y(pi/2-zoom_region[1], 0.0, Equatorial_Projection_enum)
 
  Nx = resolution
  Ny = int(<float>resolution*half_box_len_y/half_box_len_x)

  n_theta, n_phi = sample_shape(L, Method=Method)
  n_theta_north = int(n_theta/2)
  n_theta_south = n_theta-n_theta_north

  theta, phi = sample_positions(L, Grid=True, Method=Method)

  if not(rot is None):
    if len(rot) == 1:
      for i from 0 <= i < n_theta:
        for j from 0 <= j < n_phi:
          phi[i,j] += rot_angle
          if phi[i,j] > 2*pi:
            phi[i,j] -= 2*pi
          if phi[i,j] < 0:
            phi[i,j] += 2*pi

    elif len(rot) == 3:
      xx, yy, zz = s2_to_cart(theta, phi)
      xx_p, yy_p, zz_p = rot_cart_2d(xx, yy, zz, rot)
      theta, phi = cart_to_s2(xx_p, yy_p, zz_p)

  # do projection
  x_project = np.empty((theta.shape[0],theta.shape[1]), dtype=float)
  y_project = np.empty((theta.shape[0],theta.shape[1]), dtype=float)

  x_project, y_project = forward_equatorial_projection_function_array(theta, phi, Equatorial_Projection_enum)

  f_project  = np.zeros((Ny, Nx), dtype=float)
  n_points   = np.zeros((Ny, Nx), dtype=int)
  mask       = np.full((Ny, Nx), np.nan, dtype=float)

  for i in range(n_theta):
    for j in range(n_phi):

      if x_project[i,j] < half_box_len_x and x_project[i,j] > -half_box_len_x and y_project[i,j] < half_box_len_y and y_project[i,j] > -half_box_len_y:
        p = int(Nx*(x_project[i,j]+half_box_len_x)/(2.0*half_box_len_x))
        q = int(Ny*(y_project[i,j]+half_box_len_y)/(2.0*half_box_len_y))

        if np.isnan(f[i,j]):
          f_project[q,p] = np.nan
          mask[q,p]      = 0.0
        else:
          f_project[q,p] += f[i,j]
          n_points[q,p]  += 1

  if not(rot is None):
    if len(rot) == 1:
      rot_list = [rot[0],0.0,0.0]
    elif len(rot) == 3:
      rot_list = [-rot[2],-rot[1],-rot[0]]
    rot_matix = make_rotation_matrix(rot_list)

  for i in range(Nx):
    for j in range(Ny):
      if n_points[j,i] == 0:
        x_pos     = (2.*half_box_len_x*(<float>i+0.5)/<float>Nx -half_box_len_x)
        y_pos     = (2.*half_box_len_y*(<float>j+0.5)/<float>Ny -half_box_len_y)
        theta_pos = inverse_equatorial_projection_function_float_theta(x_pos, y_pos, Equatorial_Projection_enum)
        phi_pos   = inverse_equatorial_projection_function_float_phi(x_pos, y_pos, Equatorial_Projection_enum)
        # print(i, j, x_pos, y_pos, rho, theta_pos)
        if np.isnan(phi_pos):
          f_project[j,i] = np.nan
        else:
          # perform rotation
          if rot is not None:
            x_pos   = sin(theta_pos)*cos(phi_pos)
            y_pos   = sin(theta_pos)*sin(phi_pos)
            z_pos   = cos(theta_pos)
            x_p_pos = x_pos*rot_matix[0,0] + y_pos*rot_matix[0,1] + z_pos*rot_matix[0,2]
            y_p_pos = x_pos*rot_matix[1,0] + y_pos*rot_matix[1,1] + z_pos*rot_matix[1,2]
            z_p_pos = x_pos*rot_matix[2,0] + y_pos*rot_matix[2,1] + z_pos*rot_matix[2,2]

            theta_pos = atan2(sqrt(x_p_pos*x_p_pos + y_p_pos*y_p_pos),z_p_pos)
            phi_pos   = atan2(y_p_pos,x_p_pos)

          if phi_pos < 0:
            phi_pos += 2*pi
          if phi_pos > 2*pi:
            phi_pos -= 2*pi
          
          p = cy_theta_to_index(theta_pos, L, Method_enum)
          q = cy_phi_to_index(phi_pos, L, Method_enum)
          if np.isnan(f[p,q]):
            f_project[j,i] = np.nan
            mask[j,i] = 0.0
          else:
            f_project[j,i]  = f[p,q]
            n_points[j,i]  += 1
      else:
        f_project[j,i] /= <float>n_points[j,i]


  return f_project, mask

def forward_equatorial_projection_function_array(np.ndarray[ double, ndim=2, mode="c"] theta, \
  np.ndarray[ double, ndim=2, mode="c"] phi, EQUATORIAL_PROJECTION_TYPE Equatorial_Projection_enum=MERCATOR):

  cdef np.ndarray[np.float_t, ndim=2] x_project, y_project

  if Equatorial_Projection_enum == MERCATOR:
    x_project = phi
    x_project[x_project>pi] = x_project[x_project>pi] - 2*pi
    with np.errstate(divide='ignore'):
      y_project = np.log(np.tan(0.5*(pi-theta)))
    return x_project, y_project
  if Equatorial_Projection_enum == SINE:
    x_project = phi
    x_project[x_project>pi] = x_project[x_project>pi] - 2*pi
    x_project = x_project*np.sin(theta)
    y_project = pi/2-theta
    return x_project, y_project

cdef float forward_equatorial_projection_function_float_x(float theta, float phi, \
  EQUATORIAL_PROJECTION_TYPE Equatorial_Projection_enum=MERCATOR):

  cdef float x_project

  if Equatorial_Projection_enum == MERCATOR:
    x_project = phi
    if x_project > pi:
      x_project -= 2*pi
    return x_project
  if Equatorial_Projection_enum == SINE:
    x_project = phi
    if x_project > pi:
      x_project -= 2*pi
    x_project *= sin(theta)
    return x_project

cdef float forward_equatorial_projection_function_float_y(float theta, float phi,\
 EQUATORIAL_PROJECTION_TYPE Equatorial_Projection_enum=MERCATOR):

  cdef float y_project

  if Equatorial_Projection_enum == MERCATOR:
    y_project = log(tan(0.5*(pi-theta)))
    return y_project
  if Equatorial_Projection_enum == SINE:
    y_project = pi/2-theta
    return y_project

cdef float forward_equatorial_projection_function_float_x_dtheta(float theta, float phi, \
  EQUATORIAL_PROJECTION_TYPE Equatorial_Projection_enum=MERCATOR):

  cdef float dg_dtheta

  if Equatorial_Projection_enum == MERCATOR:
    return 0.0
  if Equatorial_Projection_enum == SINE:
    dg_dtheta = phi
    if dg_dtheta > pi:
      dg_dtheta -= 2*pi
    dg_dtheta *= cos(theta)
    return dg_dtheta

cdef float forward_equatorial_projection_function_float_x_dphi(float theta, float phi, \
  EQUATORIAL_PROJECTION_TYPE Equatorial_Projection_enum=MERCATOR):

  if Equatorial_Projection_enum == MERCATOR:
    return 1.0
  if Equatorial_Projection_enum == SINE:
    return sin(theta)

cdef float forward_equatorial_projection_function_float_y_dtheta(float theta, float phi,\
 EQUATORIAL_PROJECTION_TYPE Equatorial_Projection_enum=MERCATOR):

  cdef float df_dtheta

  if Equatorial_Projection_enum == MERCATOR:
    df_dtheta = -1.0/(2*sin((pi-theta)/2.0)*cos((pi-theta)/2.0))
    return df_dtheta
  if Equatorial_Projection_enum == SINE:
    df_dtheta = -1.0
    return df_dtheta

cdef float forward_equatorial_projection_function_float_y_dphi(float theta, float phi,\
 EQUATORIAL_PROJECTION_TYPE Equatorial_Projection_enum=MERCATOR):

#  cdef float df_dphi

  if Equatorial_Projection_enum == MERCATOR:
    return 0.0
  if Equatorial_Projection_enum == SINE:
    return 0.0

cdef float inverse_equatorial_projection_function_float_theta(float x, float y, \
  EQUATORIAL_PROJECTION_TYPE Equatorial_Projection_enum=MERCATOR):

  cdef float theta

  if Equatorial_Projection_enum == MERCATOR:
    theta = pi-atan(exp(y))*2.0
    return theta
  if Equatorial_Projection_enum == SINE:
    theta = pi/2-y
    return theta

cdef float inverse_equatorial_projection_function_float_phi(float x, float y, \
  EQUATORIAL_PROJECTION_TYPE Equatorial_Projection_enum=MERCATOR):

  cdef float phi

  if Equatorial_Projection_enum == MERCATOR:
    phi = x
    if phi < 0.0:
      phi += 2*pi
    return phi
  if Equatorial_Projection_enum == SINE:
    phi = x
    phi /= sin(pi/2-y)    
    if phi > pi or phi < -pi:
      return np.nan 
    if phi < 0.0:
      phi += 2*pi
    return phi

def polar_plane_to_sphere(np.ndarray[ double, ndim=2, mode="c"] image, int L, bint rot=False, Polar_Projection_enum=SP, list rotation_angles=[0.0,0.,0.], double theta_project=1.0, str Method="MW"):
  """
  Projects planar map onto a spherical MWSS map by a generalized 
  polar projection and performs random euler rotation if specified.

  Args:
    - image:
        (2D float array) image data to be augmented.
    - random_choice:
        (dictionary) dictionary of random choices.

  Returns:
    - 2D binary array of MWSS pixelized spherical image.

  Raises:
    - None
  """
  cdef METHOD_TYPE Method_enum=MW

  cdef list rotation_inverse, rotation
  cdef np.ndarray[np.float_t, ndim=2] spherical_map, spherical_count, theta_2D, phi_2D, theta_new, phi_new, x, y, z, xx, yy, zz
  cdef int i, j, theta_index, phi_index, tolerance, len_1, len_2, x_new_scal, y_new_scal, theta_max, theta_min, phi_max, phi_min
  cdef double r, phi, r_scal, theta, r_new, r_unscal, x_new, y_new, i_new, j_new
  
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

  if rot:
    rotation = rotation_angles
    rotation_inverse = []
    rotation_inverse.append(rotation[2] * -1.0)
    rotation_inverse.append(rotation[1] * -1.0)
    rotation_inverse.append(rotation[0] * -1.0)

  # spherical_map = np.zeros((L + 1, 2 * L))
  # spherical_count = np.zeros((L + 1, 2 * L))

  n_theta, n_phi = sample_shape(L, Method=Method)
  spherical_map = np.zeros((n_theta, n_phi))
  spherical_count = np.zeros((n_theta, n_phi))

  theta_2D = np.zeros((image.shape[0],image.shape[1]))
  phi_2D = np.zeros((image.shape[0],image.shape[1]))

  len_1, len_2 = image.shape[0], image.shape[1]

  r_array = np.zeros((image.shape[0], image.shape[1]))
  r2_array = np.zeros((image.shape[0], image.shape[1]))


  for i in range(image.shape[0]):
    for j in range(image.shape[1]):
        i_new = _image_index_to_unit_index(len_1, i)
        j_new = _image_index_to_unit_index(len_2, j)
        r = _planar_to_polar_coordinates_r(i_new,j_new)
        phi = _planar_to_polar_coordinates_phi(i_new,j_new)
        r_scal = r * tan(theta_project)
        theta = _polar_radius_to_theta(r_scal, Polar_Projection_enum=Polar_Projection_enum)
        theta_2D[i,j] = theta
        phi_2D[i,j] = phi
        r_array[i,j] = r_scal


  if rot:
    x, y, z = s2_to_cart(theta_2D, phi_2D)
    xx, yy, zz = rot_cart_2d(x, y, z, rotation)
    theta_2D, phi_2D = cart_to_s2(xx, yy, zz)

  # Box limits
  theta_max = 0
  theta_min = 10

  phi_max = 0
  phi_min = 10

  for i in range(image.shape[0]):
    for j in range(image.shape[1]):
        theta_index = cy_theta_to_index(theta_2D[i,j], L, Method_enum)
        phi_index = cy_phi_to_index(phi_2D[i,j], L, Method_enum)
        spherical_map[theta_index, phi_index] += image[i,j]
        spherical_count[theta_index, phi_index] += 1.0
        if theta_index > theta_max:
            theta_max = theta_index
        if theta_index < theta_min:
            theta_min = theta_index
        if phi_index > phi_max:
            phi_max = phi_index
        if phi_index < phi_min:
            phi_min = phi_index


  tolerance = int( 3 * (pi / L ) )
  theta_max += tolerance
  phi_max += tolerance
  theta_min -= tolerance
  phi_min -= tolerance

  theta_new = np.zeros((spherical_map.shape[0], spherical_map.shape[1]))
  phi_new = np.zeros((spherical_map.shape[0], spherical_map.shape[1]))

  for i in range(spherical_map.shape[0]):
    for j in range(spherical_map.shape[1]):
        if i >= theta_min and i <= theta_max and j >= phi_min and j <= phi_max:
            if spherical_count[i,j] < 1:
                theta_new[i,j] = cy_index_to_theta(i, L, Method_enum)
                phi_new[i,j] = cy_index_to_phi(j, L, Method_enum)

  if rot:
      x, y, z = s2_to_cart(theta_new, phi_new)
      xx, yy, zz = rot_cart_2d(x, y, z, rotation_inverse)
      theta_new, phi_new = cart_to_s2(xx, yy, zz)

  for i in range(spherical_map.shape[0]):
    for j in range(spherical_map.shape[1]):
        if i >= theta_min and i <= theta_max and j >= phi_min and j <= phi_max:
            if spherical_count[i,j] < 1:
                r_new = _theta_to_polar_radius(theta_new[i,j], Polar_Projection_enum=Polar_Projection_enum)
                r_unscal = r_new / tan(theta_project)
                x_new = _polar_to_planar_coordinates_x(
                    r=r_unscal, 
                    phi=phi_new[i,j]
                )
                y_new = _polar_to_planar_coordinates_y(
                    r=r_unscal, 
                    phi=phi_new[i,j]
                )
                x_new_scal = _unit_index_to_image_index(len_1, x_new)
                y_new_scal = _unit_index_to_image_index(len_2, y_new)
                if x_new_scal < image.shape[0] and x_new_scal > -1:
                    if y_new_scal < image.shape[1] and y_new_scal > -1:
                        spherical_map[i,j] = image[ x_new_scal, y_new_scal ]
            else:
                spherical_map[i,j] /= spherical_count[i,j]
        else:
            spherical_map[i,j] = 0.0


  return spherical_map

cdef inline double _image_index_to_unit_index(int len_1, int i):
  return 2.0 * (float(i)/float(len_1)  - 0.5)

cdef inline int _unit_index_to_image_index(int len_1, double i):
  return int( float(len_1) * ( float(i)/2. + 0.5) + 0.5)

cdef inline double _planar_to_polar_coordinates_r(double x, double y):
  """
  Computes polar co-ordinates from cartesian co-ordinates.

  Args:
      - x:
          (Int) cartesian x co-ordinate
      - y: 
          (Int) cartesian y co-ordinate

  Returns:
     radius r and polar angle phi as floats.

  Raises:
      - None
  """
  return sqrt( x**2 + y**2 ) 

cdef inline double _planar_to_polar_coordinates_phi(double x, double y):
  """
  Computes polar co-ordinates from cartesian co-ordinates.

  Args:
      - x:
          (Int) cartesian x co-ordinate
      - y: 
          (Int) cartesian y co-ordinate

  Returns:
     radius r and polar angle phi as floats.

  Raises:
      - None
  """
  cdef double phi
  phi = atan2(y, x)

  return phi%(2.0*pi)

cdef inline double _polar_to_planar_coordinates_x(double r, double phi):
  return r * cos(phi)

cdef inline _polar_to_planar_coordinates_y(double r, double phi):
  return r * sin(phi)

cdef inline double _polar_radius_to_theta(double r, Polar_Projection_enum=SP):
  """
  Computes the corresonding lattitude given the radius of a planar pixel
  for a variety of polar projection functions.

  Args:
      - r:    
          (Float) polar radius of planar pixel.

  Returns:
      - Theta as a float.

  Raises:
      - Exception if projection method not supported
  """

  if Polar_Projection_enum == SP:
      return 2.0 * atan(r/2.0)

  elif Polar_Projection_enum == GP:
      return atan(r)

  elif Polar_Projection_enum == OP:
      return asin(r)

  else:
      raise Exception("Projection method no supported. Try 'SP' (stereo), 'GP' (gnonmic), or 'OP' (ortho)")

cdef inline double _theta_to_polar_radius(double theta, Polar_Projection_enum=SP):
  """
  Computes the corresponding polar radius given a spherical lattitude.

  Args:
      - theta:
          (Float) Spherical lattitude in radians.

  Returns:
      - Radius as a float.

  Raises:
      - Exception if projection method not supported
  """

  if Polar_Projection_enum == SP:
      return 2.0 * tan(theta / 2.0)

  elif Polar_Projection_enum == GP:
      return tan(theta)

  elif Polar_Projection_enum == OP:
      return sin(theta)

  else:
      raise Exception("Projection method no supported. Try 'SP' (stereo), 'GP' (gnonmic), or 'OP' (ortho)")


def polar_projection_work(np.ndarray[ double, ndim=2, mode="c"] f, int L, int resolution=500, rot=None,\
                        float zoom_region=-1,  METHOD_TYPE Method_enum=MW, POLAR_PROJECTION_TYPE Polar_Projection_enum=OP):

  if rot is None or len(rot) > 1:
    dummy = 0.0
  else:
    dummy = <double> rot[0]

  if Method_enum==MW:
    Method="MW"
  elif Method_enum==MWSS:
    Method = "MWSS"
  elif Method_enum==DH:
    Method = "DH"
  elif Method_enum==GL:
    Method = "GL"

  if zoom_region < 0:
    if Polar_Projection_enum == GP:
      zoom_region = pi/4
    else:
      zoom_region = pi/2


  cdef double rot_angle = dummy

  cdef np.ndarray[np.float_t, ndim=2] theta, phi
  cdef np.ndarray[np.float_t, ndim=2] x_project, y_project
  cdef np.ndarray[np.float_t, ndim=2] xx, yy, zz, xx_p, yy_p, zz_p
  cdef np.ndarray[np.float_t, ndim=2] ortho_proj_north, mask_north, ortho_proj_south, mask_south
  cdef np.ndarray[np.float_t, ndim=2] rot_matix

  cdef list rot_list

  cdef np.ndarray[np.int_t,   ndim=2] n_points_north, n_points_south

  cdef float x_pos, y_pos, z_pos, x_p_pos, y_p_pos, z_p_pos, rho, theta_pos, phi_pos, half_box_len, max_len, tol_error
  
  cdef int n_theta, n_phi, n_theta_north, n_theta_south, i, j, p, q, i_rot, j_rot

  half_box_len = forward_projection_function_float(zoom_region, Polar_Projection_enum)
  tol_error = 1E-10

  n_theta, n_phi = sample_shape(L, Method=Method)
  n_theta_north = int(n_theta/2)
  n_theta_south = n_theta-n_theta_north

  theta, phi = sample_positions(L, Grid=True, Method=Method)

  if not(rot is None):
    if len(rot) == 1:
      for i from 0 <= i < n_theta:
        for j from 0 <= j < n_phi:
          phi[i,j] += rot_angle
          if phi[i,j] > 2*pi:
            phi[i,j] -= 2*pi
          if phi[i,j] < 0:
            phi[i,j] += 2*pi

    elif len(rot) == 3:
      xx, yy, zz = s2_to_cart(theta, phi)
      xx_p, yy_p, zz_p = rot_cart_2d(xx, yy, zz, rot)
      theta, phi = cart_to_s2(xx_p, yy_p, zz_p)

  # do projection
  x_project = np.empty((theta.shape[0],theta.shape[1]), dtype=float)
  y_project = np.empty((theta.shape[0],theta.shape[1]), dtype=float)

  x_project[theta<pi/2] = forward_projection_function_array(theta[theta<pi/2], Polar_Projection_enum)*np.cos(phi[theta<pi/2])
  y_project[theta<pi/2] = forward_projection_function_array(theta[theta<pi/2], Polar_Projection_enum)*np.sin(phi[theta<pi/2])

  x_project[theta>=pi/2] = forward_projection_function_array((pi-theta[theta>=pi/2]), Polar_Projection_enum)*np.cos(phi[theta>=pi/2])
  y_project[theta>=pi/2] = forward_projection_function_array((pi-theta[theta>=pi/2]), Polar_Projection_enum)*np.sin(phi[theta>=pi/2])

  ortho_proj_north = np.zeros((resolution, resolution), dtype=float)
  ortho_proj_south = np.zeros((resolution, resolution), dtype=float)

  n_points_north = np.zeros((resolution, resolution), dtype=int)
  n_points_south = np.zeros((resolution, resolution), dtype=int)

  mask_north   = np.empty([resolution,resolution],dtype=float)
  mask_south   = np.empty([resolution,resolution],dtype=float)

  mask_north.fill(np.nan)
  mask_south.fill(np.nan)

  for i in range(n_theta):
    for j in range(n_phi):

      if x_project[i,j] < half_box_len*(1-tol_error) and x_project[i,j] > -half_box_len*(1-tol_error) and y_project[i,j] < half_box_len*(1-tol_error) and y_project[i,j] > -half_box_len*(1-tol_error):
        p = int(resolution*(x_project[i,j]+half_box_len)/(2.0*half_box_len))
        q = int(resolution*(y_project[i,j]+half_box_len)/(2.0*half_box_len))

        if theta[i,j] < pi/2:        
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

  max_len = projection_max_len(Polar_Projection_enum)

  for i in range(resolution):
    for j in range(resolution):
      if n_points_north[i,j] == 0:
        if projection_index_to_length(i,j,resolution, half_box_len) < max_len:
          x_pos     = (2.*half_box_len*(<float>i+0.5)/<float>resolution -half_box_len)
          y_pos     = (2.*half_box_len*(<float>j+0.5)/<float>resolution -half_box_len)
          rho       = sqrt(x_pos*x_pos + y_pos*y_pos) 
          theta_pos = inverse_projection_function_float(rho, Polar_Projection_enum)
          phi_pos   = atan2(y_pos,x_pos)
          # print(i, j, x_pos, y_pos, rho, theta_pos)

          # perform rotation
          if rot is not None:
            if Polar_Projection_enum != OP:
              x_pos   = sin(theta_pos)*cos(phi_pos)
              y_pos   = sin(theta_pos)*sin(phi_pos)
            z_pos   = cos(theta_pos)
            x_p_pos = x_pos*rot_matix[0,0] + y_pos*rot_matix[0,1] + z_pos*rot_matix[0,2]
            y_p_pos = x_pos*rot_matix[1,0] + y_pos*rot_matix[1,1] + z_pos*rot_matix[1,2]
            z_p_pos = x_pos*rot_matix[2,0] + y_pos*rot_matix[2,1] + z_pos*rot_matix[2,2]

            theta_pos = atan2(sqrt(x_p_pos*x_p_pos + y_p_pos*y_p_pos),z_p_pos)
            phi_pos   = atan2(y_p_pos,x_p_pos)

          if phi_pos < 0:
            phi_pos += 2*pi
          if phi_pos > 2*pi:
            phi_pos -= 2*pi
          
          p = cy_theta_to_index(theta_pos, L, Method_enum)
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
        if projection_index_to_length(i,j,resolution, half_box_len) < max_len:
          x_pos     = (2.*half_box_len*(<float>i+0.5)/<float>resolution -half_box_len)
          y_pos     = (2.*half_box_len*(<float>j+0.5)/<float>resolution -half_box_len)
          rho       = sqrt(x_pos*x_pos + y_pos*y_pos) 
          theta_pos = pi-inverse_projection_function_float(rho, Polar_Projection_enum)
          phi_pos   = atan2(y_pos,x_pos)

          # perform rotation
          if rot is not None:
            if Polar_Projection_enum != OP:
              x_pos   = sin(theta_pos)*cos(phi_pos)
              y_pos   = sin(theta_pos)*sin(phi_pos)
            z_pos   = cos(theta_pos)
            x_p_pos = x_pos*rot_matix[0,0] + y_pos*rot_matix[0,1] + z_pos*rot_matix[0,2]
            y_p_pos = x_pos*rot_matix[1,0] + y_pos*rot_matix[1,1] + z_pos*rot_matix[1,2]
            z_p_pos = x_pos*rot_matix[2,0] + y_pos*rot_matix[2,1] + z_pos*rot_matix[2,2]

            theta_pos = atan2(sqrt(x_p_pos*x_p_pos + y_p_pos*y_p_pos),z_p_pos)
            phi_pos   = atan2(y_p_pos,x_p_pos)

          if phi_pos < 0:
            phi_pos += 2*pi
          if phi_pos > 2*pi:
            phi_pos -= 2*pi
          
          p = cy_theta_to_index(theta_pos, L, Method_enum)
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

def forward_projection_function_array(np.ndarray theta, POLAR_PROJECTION_TYPE Polar_Projection_enum):
  if Polar_Projection_enum==OP:
    return np.sin(theta)
  elif Polar_Projection_enum==GP:
    return np.tan(theta)
  elif Polar_Projection_enum==SP:
    return 2.0*np.tan(theta/2.0)

cdef inline float forward_projection_function_float(float theta, POLAR_PROJECTION_TYPE Polar_Projection_enum):
  if Polar_Projection_enum==OP:
    return sin(theta)
  elif Polar_Projection_enum==GP:
    return tan(theta)
  elif Polar_Projection_enum==SP:
    return 2.0*tan(theta/2.0)

cdef inline float forward_projection_function_prime_float(float theta, POLAR_PROJECTION_TYPE Polar_Projection_enum):
  cdef cos_theta = cos(theta)
  if Polar_Projection_enum==OP:
    cos_theta = cos(theta)
    return cos_theta
  elif Polar_Projection_enum==GP:
    cos_theta = cos(theta)
    return 1.0/(cos_theta*cos_theta)
  elif Polar_Projection_enum==SP:
    cos_theta = cos(theta/2.0)
    return 1.0/(cos_theta*cos_theta)

cdef inline float inverse_projection_function_float(float rho, POLAR_PROJECTION_TYPE Polar_Projection_enum):
  if Polar_Projection_enum==OP:
    return asin(rho)
  elif Polar_Projection_enum==GP:
    return atan(rho)
  elif Polar_Projection_enum==SP:
    return 2.0*atan(rho/2.0)

cdef inline float projection_index_to_length(int i, int j, int N, float half_box_len):
  return sqrt((2.*half_box_len*(<float>i+0.5)/<float>N -half_box_len)*(2.*half_box_len*(<float>i+0.5)/<float>N -half_box_len) \
          + (2.*half_box_len*(<float>j+0.5)/<float>N -half_box_len)*(2.*half_box_len*(<float>j+0.5)/<float>N -half_box_len))

cdef inline float projection_max_len(POLAR_PROJECTION_TYPE Polar_Projection_enum):
  if Polar_Projection_enum==OP:
    return 1.0
  elif Polar_Projection_enum==GP:
    return 1E9
  elif Polar_Projection_enum==SP:
    return 2.0


def polar_projection(f, int L, int resolution=500, rot=None,\
                        float zoom_region=-1, str Method="MW", str Projection="OP", int Spin=0):

  cdef int i, j, n_phi, n_theta
  cdef np.ndarray[np.float_t, ndim=2] f_real, f_imag
  cdef METHOD_TYPE Method_enum
  cdef POLAR_PROJECTION_TYPE Polar_Projection_enum

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

  if Projection=="OP":
    Polar_Projection_enum=OP
  elif Projection=="GP":
    Polar_Projection_enum=GP
  elif Projection=="SP":
    Polar_Projection_enum=SP
  else:
    raise ssht_input_error('Projection is not recognised, Methods are: OP, GP and SP')

  if Polar_Projection_enum==GP and zoom_region>=pi/2:
    raise ssht_input_error('zoom_region cannot be >= pi/2 for GP')
  if Polar_Projection_enum==OP and zoom_region>pi/2:
    raise ssht_input_error('zoom_region cannot be > pi/2 for OP')
  if Polar_Projection_enum==SP and zoom_region>=pi:
    raise ssht_input_error('zoom_region cannot be >= pi for SP')



  if not isinstance(f, np.ndarray):
    raise TypeError("Input not a ndarray")

  if f.ndim != 2:
    raise ssht_input_error("f must have ndim = 2")


  if f.dtype == np.float64:
    if Spin != 0:
      raise('Spin not zero but input not complex, spin signals must be complex')
    return polar_projection_work(f, L, resolution=resolution, rot=rot,\
                        zoom_region=zoom_region, Method_enum=Method_enum, Polar_Projection_enum=Polar_Projection_enum)
  elif f.dtype == complex:
    n_theta, n_phi = sample_shape(L,Method=Method)
    f_real = np.empty((n_theta,n_phi), dtype=np.float_)
    f_imag = np.empty((n_theta,n_phi), dtype=np.float_)
    for i in range(n_theta):
      for j in range(n_phi):
        f_real[i,j] = f[i,j].real
        f_imag[i,j] = f[i,j].imag

    north_plot_real, mask_north_real, south_plot_real, mask_south_real\
                         = polar_projection_work(f_real, L, resolution=resolution, rot=rot,\
                        zoom_region=zoom_region, Method_enum=Method_enum, Polar_Projection_enum=Polar_Projection_enum)
    north_plot_imag, mask_north_imag, south_plot_imag, mask_south_imag\
                        = polar_projection_work(f_imag, L, resolution=resolution, rot=rot,\
                        zoom_region=zoom_region, Method_enum=Method_enum, Polar_Projection_enum=Polar_Projection_enum)
    if Spin != 0:
      rotation_angle_north, rotation_angle_south = polar_projection_rotation_array(resolution, \
                          Projection=Projection, rot=rot, zoom_region=zoom_region)
      north_plot_real_keep = north_plot_real
      south_plot_real_keep = south_plot_real
      north_plot_real = north_plot_real*np.cos(Spin*rotation_angle_north) - north_plot_imag*np.sin(Spin*rotation_angle_north)
      north_plot_imag = north_plot_imag*np.cos(Spin*rotation_angle_north) + north_plot_real_keep*np.sin(Spin*rotation_angle_north)

      south_plot_real = south_plot_real*np.cos(Spin*rotation_angle_south) - south_plot_imag*np.sin(Spin*rotation_angle_south)
      south_plot_imag = south_plot_imag*np.cos(Spin*rotation_angle_south) + south_plot_real_keep*np.sin(Spin*rotation_angle_south)
    return north_plot_real, mask_north_real, south_plot_real, mask_south_real,\
           north_plot_imag, mask_north_imag, south_plot_imag, mask_south_imag
  else:
    raise ssht_input_error("f dtype must be float or complex")

def polar_projection_rotation_array(int resolution, str Projection="OP", rot=None, float zoom_region=-1):

  cdef np.ndarray[np.float_t, ndim=2] rotation_angle_north, rotation_angle_south
  cdef POLAR_PROJECTION_TYPE Polar_Projection_enum

  if Projection=="OP":
    Polar_Projection_enum=OP
  elif Projection=="GP":
    Polar_Projection_enum=GP
  elif Projection=="SP":
    Polar_Projection_enum=SP
  else:
    raise ssht_input_error('Projection is not recognised, Methods are: OP, GP and SP')

  if zoom_region < 0:
    if Polar_Projection_enum == GP:
      zoom_region = pi/4
    else:
      zoom_region = pi/2


  cdef np.ndarray[np.float_t, ndim=2] rot_matix

  cdef list rot_list

  cdef double x_pos_plane, y_pos_plane, z_pos_plane
  cdef double rho, theta_pos, phi_pos, half_box_len, max_len
  cdef double delta_x, delta_y, delta_z, delta_rho, delta_phi, delta_x_plane, delta_y_plane, r_xy
  
  cdef int n_theta, n_phi, n_theta_north, n_theta_south, i, j, p, q, i_rot, j_rot

  rotation_angle_north = np.full((resolution,resolution),np.nan, dtype=float)
  rotation_angle_south = np.full((resolution,resolution),np.nan, dtype=float)

  half_box_len = forward_projection_function_float(zoom_region, Polar_Projection_enum)
 
  if not(rot is None):
    if len(rot) == 1:
      rot_list = [-rot[0],0.0,0.0]
    elif len(rot) == 3:
      rot_list = rot
    rot_matix = make_rotation_matrix(rot_list)
    delta_x = rot_matix[0,2]
    delta_y = rot_matix[1,2]
    delta_z = rot_matix[2,2]
  else:
    delta_x = 0.0
    delta_y = 0.0
    delta_z = 1.0

  max_len = projection_max_len(Polar_Projection_enum)

  for i in range(resolution):
    for j in range(resolution):
        if projection_index_to_length(i,j,resolution, half_box_len) < max_len:
          x_pos_plane     = (2.*half_box_len*(<float>i+0.5)/<float>resolution -half_box_len)
          y_pos_plane     = (2.*half_box_len*(<float>j+0.5)/<float>resolution -half_box_len)
          rho       = sqrt(x_pos_plane*x_pos_plane + y_pos_plane*y_pos_plane) 
          theta_pos = inverse_projection_function_float(rho, Polar_Projection_enum)
          phi_pos   = atan2(y_pos_plane,x_pos_plane)

          delta_rho  = sin(phi_pos)*delta_y+cos(phi_pos)*delta_x - tan(theta_pos)*delta_z
          delta_rho *= cos(theta_pos)
          delta_rho *= forward_projection_function_prime_float(theta_pos, Polar_Projection_enum)

          delta_phi  = -tan(phi_pos)*delta_x + delta_y
          delta_phi *= cos(phi_pos)/sin(theta_pos)

          delta_x_plane = cos(phi_pos)*delta_rho - rho*sin(phi_pos)*delta_phi
          delta_y_plane = sin(phi_pos)*delta_rho + rho*cos(phi_pos)*delta_phi

          # calculate angle
          rotation_angle_north[i,j] = -atan2(delta_x_plane, delta_y_plane)

        if projection_index_to_length(i,j,resolution, half_box_len) < max_len:
          x_pos_plane     = (2.*half_box_len*(<float>i+0.5)/<float>resolution -half_box_len)
          y_pos_plane     = (2.*half_box_len*(<float>j+0.5)/<float>resolution -half_box_len)
          rho       = sqrt(x_pos_plane*x_pos_plane + y_pos_plane*y_pos_plane) 
          theta_pos = pi-inverse_projection_function_float(rho, Polar_Projection_enum)
          phi_pos   = atan2(y_pos_plane,x_pos_plane)

          delta_rho  = sin(phi_pos)*delta_y+cos(phi_pos)*delta_x - tan(theta_pos)*delta_z
          delta_rho *= cos(theta_pos)
          delta_rho *= -forward_projection_function_prime_float(pi-theta_pos, Polar_Projection_enum)

          delta_phi  = -tan(phi_pos)*delta_x + delta_y
          delta_phi *= cos(phi_pos)/sin(theta_pos)

          delta_x_plane = cos(phi_pos)*delta_rho - rho*sin(phi_pos)*delta_phi
          delta_y_plane = sin(phi_pos)*delta_rho + rho*cos(phi_pos)*delta_phi

          # calculate angle
          rotation_angle_south[i,j] = atan2(delta_x_plane, delta_y_plane)



  return rotation_angle_north, rotation_angle_south



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

  MAX_ITERATIONS = 1e5
  TOL = 1e-10

  #% Convert theta to longitude.
  thetas, phis = theta_phi_to_ra_dec(thetas,phis)

  t = thetas
  cdef bint inaccurate = True
  cdef int count = 0

  while(inaccurate):
    count += 1
    dt = (t + np.sin(t) - pi*np.sin(thetas)) / (1 + np.cos(t))
    t = t - dt
    if(np.max(np.abs(dt)) < TOL or count > MAX_ITERATIONS):
      inaccurate = False

  t = t/2
  x = 2 * np.sqrt(2) / pi * phis * np.cos(t)
  y = np.sqrt(2) * np.sin(t)

  return (x, y)

#----------------------------------------------------------------------------------------------------#

# dl functions

def dl_beta_recurse(np.ndarray[ double, ndim=2, mode="c"] dl not None, double beta, int L,\
            int el, np.ndarray[ double, ndim=1, mode="c"] sqrt_tbl not None,\
            np.ndarray[ double, ndim=1, mode="c"] signs not None):

  cdef ssht_dl_size_t dl_size=SSHT_DL_HALF
  
  ssht_dl_beta_risbo_half_table(<double*> np.PyArray_DATA(dl), beta, L, dl_size,\
                                el, <double*> np.PyArray_DATA(sqrt_tbl),<double*> np.PyArray_DATA(signs))

  return dl

def dln_beta_recurse(np.ndarray[ double, ndim=1, mode="c"] dl not None,\
            np.ndarray[ double, ndim=1, mode="c"] dlm1 not None, double beta,\
            int L, int el, int n, np.ndarray[ double, ndim=1, mode="c"] sqrt_tbl not None,\
            np.ndarray[ double, ndim=1, mode="c"] signs not None):

  ssht_dl_beta_kostelec_halfline_table(<double*> np.PyArray_DATA(dlm1),\
                                        <double*> np.PyArray_DATA(dl),\
                                        beta, L, n, el,\
                                        <double*> np.PyArray_DATA(sqrt_tbl),\
                                        <double*> np.PyArray_DATA(signs))

  return dlm1

def generate_dl(double beta, int L):

  cdef np.ndarray[np.float_t, ndim=3] dl_array = np.zeros((L, 2*L-1, 2*L-1), dtype=np.float_)
  cdef np.ndarray[np.float_t, ndim=2] dl_dummy = np.zeros((2*L-1, 2*L-1), dtype=np.float_)

  cdef np.ndarray[np.float_t, ndim=1] sqrt_tbl = np.sqrt(np.arange(0,2*(L-1)+1, dtype=np.float_))
  cdef np.ndarray[np.float_t, ndim=2] signs = np.ones((L+1,1), dtype=np.float_)

  cdef int i, j, el, offset_m=L-1

  for i in range(1,L+1,2):
    signs[i] = -1 

  
  cdef ssht_dl_size_t dl_size=SSHT_DL_HALF
  
  # do recursion
  #el = 0 first
  ssht_dl_beta_risbo_half_table(<double*> np.PyArray_DATA(dl_dummy), beta, L, dl_size,\
                                           0, <double*> np.PyArray_DATA(sqrt_tbl), <double*> np.PyArray_DATA(signs))

  dl_array[0,offset_m,offset_m] = dl_dummy[offset_m,offset_m]

  for el in range(1,L):
    ssht_dl_beta_risbo_half_table(<double*> np.PyArray_DATA(dl_dummy), beta, L, dl_size,\
                                           el, <double*> np.PyArray_DATA(sqrt_tbl), <double*> np.PyArray_DATA(signs))
    for i in range(-el,el+1):
      # print(i)
      for j in range(-el,el+1):
        dl_array[el,offset_m+i,offset_m+j] = dl_dummy[offset_m+i,offset_m+j]

  return dl_array

def generate_dl_Mmax(double beta, int L, int M):

  cdef np.ndarray[np.float_t, ndim=3] dl_array = np.zeros((L, 2*M-1, 2*M-1), dtype=np.float_)
  cdef np.ndarray[np.float_t, ndim=2] dl_dummy = np.zeros((2*L-1, 2*L-1), dtype=np.float_)

  cdef np.ndarray[np.float_t, ndim=1] sqrt_tbl = np.sqrt(np.arange(0,2*(L-1)+1, dtype=np.float_))
  cdef np.ndarray[np.float_t, ndim=2] signs = np.ones((L+1,1), dtype=np.float_)

  cdef int i, j, el, offset_m=M-1, offset_m_dummy=L-1

  for i in range(1,L+1,2):
    signs[i] = -1 

  
  cdef ssht_dl_size_t dl_size=SSHT_DL_HALF
  
  # do recursion
  #el = 0 first
  ssht_dl_beta_risbo_half_table(<double*> np.PyArray_DATA(dl_dummy), beta, L, dl_size,\
                                           0, <double*> np.PyArray_DATA(sqrt_tbl), <double*> np.PyArray_DATA(signs))

  dl_array[0,offset_m,offset_m] = dl_dummy[offset_m_dummy,offset_m_dummy]

  for el in range(1,L):
    # print(el)
    ssht_dl_beta_risbo_half_table(<double*> np.PyArray_DATA(dl_dummy), beta, L, dl_size,\
                                           el, <double*> np.PyArray_DATA(sqrt_tbl), <double*> np.PyArray_DATA(signs))
    for i in range(-M+1,M):
      # print(i)
      for j in range(-M+1,M):
        dl_array[el,offset_m+i,offset_m+j] = dl_dummy[offset_m_dummy+i,offset_m_dummy+j]

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


def rotate_flms(
  np.ndarray[ double complex, ndim=1, mode="c"] f_lm not None,
  double alpha,
  double beta,
  double gamma,
  int L,
  dl_array_in=None,
  M_in=None,
  bint Axisymmetric=False,
  bint Keep_dl=False,
  backend: str = "SSHT",
  **kwargs,
):
  from pyssht.ducc_interface import rotate_flms as _ducc0_rotate_flms
  from pyssht.parameters import method, Ducc

  params = method(backend=backend, **kwargs)
  if isinstance(params, Ducc):
     return _ducc0_rotate_flms(f_lm, alpha, beta, gamma, L, nthreads=params.nthreads)

  cdef np.ndarray[np.float_t, ndim=3] dl_array
  cdef np.ndarray[complex, ndim=1] alpha_array, gamma_array 
  cdef np.ndarray[complex, ndim=1] f_lm_rotated
  cdef int M, el, m

  cdef int index = 0, ind, n_max, n
  cdef complex Dlmn=0

  if dl_array_in is None:
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
                   * <complex> gamma_array[n+L-1]
            if Axisymmetric:
                ind = el
            else:
                ind = cy_elm2ind(el,n)

            f_lm_rotated[index] = <double complex> f_lm_rotated[index] + \
                <complex> Dlmn * <complex> f_lm[ind]

        index = index + 1

  if Keep_dl:
    return f_lm_rotated, dl_array
  else:
    return f_lm_rotated


#----------------------------------------------------------------------------------------------------#

# other usefull functions

def gaussian_smoothing(np.ndarray[ double complex, ndim=1, mode="c"] f_lm not None, int L, sigma_in=None, bl_in = None):

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
      bl[el] = exp(-<float>el*<float>el*sigma*sigma)

  if bl_in is not None:
    bl = bl_in

  for el in range(L):
    for m in range(-el,el+1):
      index = cy_elm2ind(el, m)
      fs_lm[index] = f_lm[index]*bl[el]

  return fs_lm

def kostelec_ylm(thetas, phis, int L, int Spin=0):
  '''
  compute the ylm with Kostelec recursion
  '''
  # initialise
  ylm = np.zeros((L * L, thetas.size), dtype=complex)

  # Precompute data.
  sqrt_tbl = np.sqrt(range(2 * L))
  signs = np.ones(L+1)
  signs[1::2] = -1

  # Compute spherical harmonics from Wigner functions.
  for itheta, theta in enumerate(thetas):
    dln = np.zeros(L)
    dlnm1 = np.zeros(L)
    for el in range(abs(Spin), L):
      tmp = dln
      dln = dln_beta_recurse(dln, dlnm1, theta, L, el, -Spin, sqrt_tbl, signs)
      dlnm1 = tmp
      ######## m = 0 ########
      ind = cy_elm2ind(el, 0)
      ylm[ind][itheta] = signs[abs(Spin)] * np.sqrt((2*el+1)/(4*np.pi)) * dln[0]
      for m in range(1, el + 1):
        ind_pm = cy_elm2ind(el, m)
        ylm[ind_pm][itheta] = signs[abs(Spin)] * np.sqrt((2*el+1)/(4*np.pi))\
        * dln[m] * np.exp(1j * m * phis[itheta])
        ind_nm = cy_elm2ind(el, -m)
        ylm[ind_nm][itheta] = (-1) ** m * np.conj(ylm[ind_pm][itheta])

  return ylm

def risbo_ylm(thetas, phis, int L, int Spin=0):
  '''
  compute the ylm with Risbo recursion
  '''
  # initialise
  ylm = np.zeros((L * L, thetas.size), dtype=complex)

  # Precompute data.
  sqrt_tbl = np.sqrt(range(2 * L))
  signs = np.ones(L+1)
  signs[1::2] = -1

  # Compute spherical harmonics from Wigner functions.
  for itheta, theta in enumerate(thetas):
    dl = np.zeros((2 * L - 1, 2 * L - 1))
    for el in range(abs(Spin), L):
      dl = dl_beta_recurse(dl, theta, L, el, sqrt_tbl, signs)
      ######## m = 0 ########
      ind = cy_elm2ind(el, 0)
      ylm[ind][itheta] = signs[abs(Spin)] * np.sqrt((2*el+1)/(4*np.pi)) * dl[L - 1][-Spin + L - 1]
      for m in range(1, el + 1):
        ind_pm = cy_elm2ind(el, m)
        ylm[ind_pm][itheta] = signs[abs(Spin)] * np.sqrt((2*el+1)/(4*np.pi))\
        * dl[m + L - 1][-Spin + L - 1] * np.exp(1j * m * phis[itheta])
        ind_nm = cy_elm2ind(el, -m)
        ylm[ind_nm][itheta] = (-1) ** m * np.conj(ylm[ind_pm][itheta])

  return ylm

def numerical_recipes_ylm(thetas, phis, int L, int Spin=0):
  '''
  compute the ylm with numerical recipes recursion
  '''
  if Spin != 0:
      raise ssht_input_error('Non-zero spin not supported for NumericalRecipes ylm recursion.')

  # initialise
  ylm = np.zeros((L * L, thetas.size), dtype=complex)

  for itheta, theta in enumerate(thetas):
    x = np.cos(theta)

    ######## m = 0 ########
    # Compute Pmm.
    pmm = 1
    somx2 = np.sqrt((1 - x) * (1 + x))
    fact = 1
    ind = cy_elm2ind(0, 0)
    c1 = 1 / (4 * np.pi)
    C = np.sqrt(c1)
    ylm[ind][itheta]  = C * pmm

    # Compute Pm,m+1.
    if L - 1 != 0:
      pmmp1 = x * pmm
      ind = cy_elm2ind(1, 0)
      c1 = 3 / (4 * np.pi)
      C = np.sqrt(c1)
      ylm[ind][itheta] = C * pmmp1

    # Compute Pm,l for l > m+1.
    for el in range(2, L):
      pll = (x * (2 * el - 1) * pmmp1 - (el - 1) * pmm) / el
      pmm = pmmp1
      pmmp1 = pll
      ind = cy_elm2ind(el, 0)
      c1 = (2 * el + 1) / (4 * np.pi)
      C = np.sqrt(c1)
      ylm[ind][itheta] = C * pll

    ######## m != 0 ########
    for m in range(1, L):
      # Compute Pmm.
      pmm = 1
      somx2 = np.sqrt((1 - x) * (1 + x))
      fact = 1
      for i in range(1, m + 1):
        pmm = -pmm * fact * somx2
        fact = fact + 2.0
      ind_pm = cy_elm2ind(m, m)
      c1 = (2 * m + 1) / (4 * np.pi)
      C = np.sqrt(c1 * factorial(
        m - abs(m), exact=False) / factorial(
          m + abs(m), exact=False))
      ylm[ind_pm][itheta] = C * pmm * np.exp(1j * m * phis[itheta])
      ind_nm = cy_elm2ind(m, -m)
      ylm[ind_nm][itheta] = (-1) ** m * np.conj(ylm[ind_pm][itheta])

      # Compute Pm,m+1.
      if m != L - 1:
        pmmp1 = x * (2 * m + 1) * pmm
        ind_pm = cy_elm2ind(m + 1, m)
        c1 = (2 * (m + 1) + 1) / (4 * np.pi)
        C = np.sqrt(c1 * factorial(
          m + 1 - abs(m), exact=False) / factorial(
            m + 1 +abs(m), exact=False))
        ylm[ind_pm][itheta] = C * pmmp1 * np.exp(1j * m * phis[itheta])
        ind_nm = cy_elm2ind(m + 1, -m)
        ylm[ind_nm][itheta] = (-1) ** m * np.conj(ylm[ind_pm][itheta])

      # Compute Pm,l for l > m+1.
      for el in range(m + 2, L):
        pll = (x * (2 * el - 1) * pmmp1 - (el + m - 1) * pmm) / (el - m)
        pmm = pmmp1
        pmmp1 = pll
        ind_pm = cy_elm2ind(el, m)
        c1 = (2 * el + 1) / (4 * np.pi)
        C = np.sqrt(c1 * factorial(
          el - abs(m), exact=False) / factorial(
            el + abs(m), exact=False))
        ylm[ind_pm][itheta] = C * pll * np.exp(1j * m * phis[itheta])
        ind_nm = cy_elm2ind(el, -m)
        ylm[ind_nm][itheta] = (-1) ** m * np.conj(ylm[ind_pm][itheta])

  return ylm

def create_ylm(thetas, phis, int L, int Spin=0, str recursion='Risbo'):
  # check if thetas is number or array and reshape if not
  if isinstance(thetas, (int, float)):
    theta_m, theta_n = 1, 1
    thetas = np.array([thetas])
  elif thetas.ndim == 1:
    theta_m, theta_n = len(thetas), 1
  else:
    theta_m, theta_n = np.shape(thetas)
    thetas = thetas.flatten()

  # check if phis is number or array and reshape if not
  if isinstance(phis, (int, float)):
    phi_m, phi_n = 1, 1
    phis = np.array([phis])
  elif phis.ndim == 1:
    phi_m, phi_n = len(phis), 1
  else:
    phi_m, phi_n = np.shape(phis)
    phis = phis.flatten()

  # Check size theta and phi identical.
  if theta_m != phi_m or theta_n != phi_n:
    raise ssht_input_error('Inconsistent theta and phi data.')

  if recursion == 'Kostelec':
    ylm = kostelec_ylm(thetas, phis, L, Spin)
  elif recursion == 'Risbo':
    ylm = risbo_ylm(thetas, phis, L, Spin)
  elif recursion == 'NumericalRecipes':
    ylm = numerical_recipes_ylm(thetas, phis, L, Spin)
  else:
    raise ssht_input_error('Invalid recursion method.')

  # Reshape output data.
  ylm = ylm.reshape(-1, theta_m, theta_n)

  return ylm
