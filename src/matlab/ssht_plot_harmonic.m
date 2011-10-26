function ssht_plot_harmonic(flm, L)
% ssht_plot_harmonic - Plot spherical harmonic coefficients
%
% Plots spherical harmonic coefficients.
%
% Default usage is given by
%
%   ssht_plot_harmonic(flm, L)
%
% where flm is the spherical harmonic coefficients of the function and L is
% the harmonic band-limit.
%
% Author: Jason McEwen (www.jasonmcewen.org)

% SSHT package to perform spin spherical harmonic transforms
% Copyright (C) 2011  Jason McEwen
% See LICENSE.txt for license details

flm_plot = zeros(2*L-1, L);

% Impose reality on flms.
for el = 0:L-1   
   for m = -el:el
      ind = ssht_elm2ind(el, m);
      flm_plot(m+L+1, el+1) = flm(ind);                  
   end
end

imagesc(flm_plot);
