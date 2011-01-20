
#ifndef SSHT_DL
#define SSHT_DL


/*!
 * Size of dl plane required.
 */
typedef enum {SSHT_DL_QUARTER = 0, SSHT_DL_HALF} ssht_dl_size_t;



double* ssht_dl_calloc(int L, ssht_dl_size_t dl_size);
int ssht_dl_get_mmoffset(int L, ssht_dl_size_t dl_size);
int ssht_dl_get_mmstride(int L, ssht_dl_size_t dl_size);
void ssht_dl_halfpi_trapani_eighth_table(double *dl, int L, 
					 ssht_dl_size_t dl_size,
					 int el, double *sqrt_tbl);


#endif
