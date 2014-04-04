function ssht_idl_test

  soname = '/Users/bl/Downloads/ssht/lib/c/libssht.dylib'

  dl = [0.5, 0.2]
  theta = 0.5
  L = 10
  el = 3
  sqrt_tble = 1.0 + dblarr(2*L+1)
  signs = 1.0 + dblarr(L+1)
  r = call_external(soname, 'ssht_idl_dl', dl, theta, L, el, sqrt_tble, signs, /CDECL)
  print, 'Return value :', r
  return, dl
end
