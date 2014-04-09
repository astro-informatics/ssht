function ssht_idl_test

  soname = '../../lib/c/libssht.dylib'
  
  big_l = 10L ; l_max + 1
  dl = DBLARR((2L * big_l - 1L) ^ 2)
  theta = 0.5D
  sqrt_tbl = SQRT(DINDGEN(2L * (big_l - 1L) + 2L))
  signs = DBLARR(big_l + 1L)
  FOR em = 0, big_L - 1L, 2 DO BEGIN
     
     signs(em) = 1.0D
     signs(em + 1) = -1.0D
     
  ENDFOR
  FOR el = 0, big_l - 1L DO BEGIN
     
     PRINT, el
     r = CALL_EXTERNAL(soname, 'ssht_idl_dl', dl, theta, big_l, el, sqrt_tbl, signs, /CDECL)
  
  ENDFOR
  print, 'Return value :', r
  return, dl
  
end
