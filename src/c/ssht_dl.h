
#ifndef SSHT_DL
#define SSHT_DL


/*! Size of dl plane required (for memory access). */
typedef enum {SSHT_DL_QUARTER = 0, SSHT_DL_HALF, SSHT_DL_FULL} ssht_dl_size_t;


double* ssht_dl_calloc(int L, ssht_dl_size_t dl_size);
int ssht_dl_get_offset(int L, ssht_dl_size_t dl_size);
int ssht_dl_get_stride(int L, ssht_dl_size_t dl_size);
void ssht_dl_beta_risbo_full_table(double *dl, double beta, int L, 
				   ssht_dl_size_t dl_size,
				   int el, double *sqrt_tbl);
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

#endif
