% ssht_demo1 - Run demo1
%
% Simple demo to compute inverse and forward transform of complex scalar
% function, using simplest interface with default options. 
%
% Default usage is given by
%
%   ssht_demo1
%
% Author: Jason McEwen (www.jasonmcewen.org)

clear all;

% Define parameters.
L = 64

% Generate random flms (of complex signal).
flm = zeros(L^2,1);
flm = rand(size(flm)) + sqrt(-1)*rand(size(flm));
flm = 2.*(flm - (1+sqrt(-1))./2);

% Compute inverse then forward transform.
f = ssht_inverse(flm, L);
flm_syn = ssht_forward(f, L);

% Compute max error in harmonic space.
maxerr = max(abs(flm_syn - flm))

% Compute sampling grids.
[thetas, phis, n, ntheta, nphi] = ssht_sampling(L, 'Grid', true);

% Plot function on sphere.
figure;
ssht_plot_sphere(abs(f), L);


