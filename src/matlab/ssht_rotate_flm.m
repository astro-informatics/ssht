function [flm_rotated] = ssht_rotate_flm(flm, alpha, beta, gamma)
% ssht_rotate_flm - Rotate a function in harmonic space
%
% Rotates a function on the sphere in harmonic space, where
% the rotation is specified by a point (alpha, beta, gamma)
% in SO(3).
%
% Default usage is given by
%
%   [flm_rotated] = ssht_rotate_flm(flm, alpha, beta, gamma)
%
% where flm are the harmonic coefficients of the function to be
% rotated, alpha, beta and gamma are the Euler angles specifying the
% rotation, and flm_rotated are the harmonic coefficients of the resulting
% function after rotation.  The band-limit L is implied by the size of the 
% flm, which should be L*L.
%
% Authors: Jason McEwen (www.jasonmcewen.org)

% SSHT package to perform spin spherical harmonic transforms
% Copyright (C) 2011-2017 Jason McEwen
% See LICENSE.txt for license details

% Parse arguments.
p = inputParser;
p.addRequired('flm', @isnumeric);          
p.addRequired('alpha', @isnumeric);          
p.addRequired('beta', @isnumeric);          
p.addRequired('gamma', @isnumeric);          

L = sqrt(length(flm));
 
% Rotate harmonic coefficients.
flm_rotated = zeros(L^2, 1);
dl = zeros(2*L-1,2*L-1);
for el = 0:L-1
  dl = ssht_dl(dl, L, el, beta);
  
  for m = -el:el
    ind_m = ssht_elm2ind(el, m);
    
    for n = -el:el
      
      ind_n = ssht_elm2ind(el, n);
      flm_rotated(ind_m) = flm_rotated(ind_m) ...
        + exp(-1i .* m .* alpha) .* dl(m+L, n+L) ...
        .* exp(-1i .* n .* gamma) .* flm(ind_n);
      
    end
    
  end
end
