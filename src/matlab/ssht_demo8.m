% ssht_demo8 - Run demo8
%
% Plot spherical harmonic functions on the sphere.
% 
% Default usage is given by
%
%   ssht_demo8
%
% Author: Jason McEwen (www.jasonmcewen.org)

% SSHT package to perform spin spherical harmonic transforms
% Copyright (C) 2011  Jason McEwen
% See LICENSE.txt for license details

clear all;
close all;

% Define parameters.
L = 64;
el = 4
m = 3
type = 'colour';

% Generate spherical harmonics.
flm = zeros(L^2,1);
ind = ssht_elm2ind(el, m);
flm(ind) = 1.0;

% Compute function on the sphere.
f = ssht_inverse(complex(real(flm), imag(flm)), L);

% Compute sampling grids.
[thetas, phis, n, ntheta, nphi] = ssht_sampling(L, 'Grid', true);

% Plot function on sphere.
figure(1);
subplot(2,1,1)
ssht_plot_sphere(real(f), L, 'Type', type, ...
    'ColourBar', true, 'Lighting', true);
subplot(2,1,2)
ssht_plot_sphere(imag(f), L, 'Type', type, ...
    'ColourBar', true, 'Lighting', true);

