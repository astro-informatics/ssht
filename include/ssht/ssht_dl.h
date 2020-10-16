// SSHT package to perform spin spherical harmonic transforms
// Copyright (C) 2011  Jason McEwen
// See LICENSE.txt for license details


#ifndef SSHT_DL
#define SSHT_DL


#ifdef __cplusplus
extern "C"{
#endif

/*! Size of dl plane required (for memory access). */
typedef enum {SSHT_DL_QUARTER = 0, 
	      SSHT_DL_QUARTER_EXTENDED, 
	      SSHT_DL_HALF, 
	      SSHT_DL_FULL} ssht_dl_size_t;

/*! Recursion to use to compute dl plane for MW methods. */
typedef enum {SSHT_DL_RISBO = 0, 
	      SSHT_DL_TRAPANI} ssht_dl_method_t;


double* ssht_dl_calloc(int L, ssht_dl_size_t dl_size);
int ssht_dl_get_offset(int L, ssht_dl_size_t dl_size);
int ssht_dl_get_stride(int L, ssht_dl_size_t dl_size);

void ssht_dl_beta_risbo_full_table(double *dl, double beta, int L, 
				   ssht_dl_size_t dl_size,
				   int el, double *sqrt_tbl);
void ssht_dl_beta_risbo_half_table(double *dl, double beta, int L, 
					 				   ssht_dl_size_t dl_size,
					 				   int el, double *sqrt_tbl,
					 				   double *signs);
void ssht_dl_beta_risbo_eighth_table(double *dl, double beta, int L, 
				     ssht_dl_size_t dl_size,
				      int el, double *sqrt_tbl,
				      double *signs);
void ssht_dl_beta_risbo_fill_eighth2quarter_table(double *dl, 
						  double *dl8,
						  int L,
						  ssht_dl_size_t dl_size,
						  ssht_dl_size_t dl8_size,
						  int el, 
						  double *signs);

void ssht_dl_beta_kostelec_full_table(double *dlm1p1, double *dl, 
				 double beta, int L, 
				 ssht_dl_size_t dl_size,
				 int el, 
				 double *sqrt_tbl, double *signs);
void ssht_dl_beta_kostelec_line_table(double *dlm1p1_line, double *dl_line, 
				      double beta, int L, int mm, int el, 
				      double *sqrt_tbl, double *signs);
void ssht_dl_beta_kostelec_halfline_table(double *dlm1p1_line, double *dl_line, 
					  double beta, int L, int mm, int el, 
					  double *sqrt_tbl, double *signs);

void ssht_dl_halfpi_trapani_eighth_table(double *dl, int L, 
					 ssht_dl_size_t dl_size,
					 int el, double *sqrt_tbl);
void ssht_dl_halfpi_trapani_quarter_table(double *dl, int L, 
					  ssht_dl_size_t dl_size,
					  int el, double *sqrt_tbl);
void ssht_dl_halfpi_trapani_fill_eighth2righthalf_table(double *dl, int L,
							ssht_dl_size_t dl_size,
							int el, double *signs);
void ssht_dl_halfpi_trapani_fill_eighth2quarter_table(double *dl, int L,
						     ssht_dl_size_t dl_size,
						     int el, double *signs);

#ifdef __cplusplus
}
#endif

#endif
