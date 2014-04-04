
#include "ssht_dl.h"

/*
 * IDL integer types. For historical reasons, we use UCHAR for TYP_BYTE
 * instead of defining an IDL_BYTE type.
 */
#if defined(ALPHA_OSF) || defined(SUN_64) || defined(LINUX_X86_64) || defined(HPUX_64) || defined(IRIX_64) || defined(AIX_64)
#define IDL_SIZEOF_C_LONG 8
#else
#define IDL_SIZEOF_C_LONG 4
#endif
#if (IDL_SIZEOF_C_LONG == 8) || defined(MSWIN_64)
#define IDL_SIZEOF_C_PTR  8
#else
#define IDL_SIZEOF_C_PTR  4
#endif
typedef short IDL_INT;
typedef unsigned short IDL_UINT;
#if IDL_SIZEOF_C_LONG == 8
typedef int IDL_LONG;
typedef unsigned int IDL_ULONG;
#elif IDL_SIZEOF_C_LONG == 4
typedef long IDL_LONG;
typedef unsigned long IDL_ULONG;
#else
#error "IDL_LONG not defined --- unexpected value of IDL_SIZEOF_C_LONG"
#endif


int ssht_idl_dl(int argc, void* argv[])  
{  
  if(argc != 6) return 0;  
  double *dl = (double *) argv[0];
  double *theta = (double *) argv[1];
  IDL_INT *L = (IDL_INT *) argv[2];
  IDL_INT *el = (IDL_INT *) argv[3];
  double *sqrt_tbl = (double *) argv[4];
  double *signs = (double *) argv[5];

  ssht_dl_beta_risbo_half_table(dl, *theta, *L, SSHT_DL_FULL, 
		*el, sqrt_tbl, signs);

  return 1;  
}    

