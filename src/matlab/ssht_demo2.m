% ssht_demo2 - Run demo2
%
% Simple demo to compute inverse and forward transform of real scalar
% function, using simplest interface with default options. 
%
% Default usage is given by
%
%   ssht_demo2
%
% Author: Jason McEwen (www.jasonmcewen.org)

clear all;

% Define parameters.
L = 64

% Generate random flms (of complex signal).
flm = zeros(L^2,1);
flm = rand(size(flm)) + sqrt(-1)*rand(size(flm));
flm = 2.*(flm - (1+sqrt(-1))./2);

% Impose reality on flms.
for el = 0:L-1
   m = 0;
   ind = ssht_elm2ind(el, m);
   flm(ind) = real(flm(ind));
   for m = 1:el
      ind_pm = ssht_elm2ind(el, m);
      ind_nm = ssht_elm2ind(el, -m);
      flm(ind_nm) = (-1)^m * conj(flm(ind_pm));
   end
end

% Compute inverse then forward transform.
f = ssht_inverse(flm, L, 'Reality', true);
flm_syn = ssht_forward(f, L, 'Reality', true);

% Compute max error in harmonic space.
maxerr = max(abs(flm_syn - flm))

% Compute sampling grids.
[thetas, phis, n, ntheta, nphi] = ssht_sampling(L, 'Grid', true);

% Plot function on sphere.
figure;
ssht_plot_sphere(f, L);


