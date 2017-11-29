% ssht_demo_rotation - Run rotation demo
%
% Plot spherical harmonic function on the sphere, then
% rotate it and plot the rotation, too.
%
% Default usage is given by
%
%   ssht_demo_rotation
%
% Authors: Martin Büttner (m.buettner.d@gmail.com)
%          Jason McEwen (www.jasonmcewen.org)

% SSHT package to perform spin spherical harmonic transforms
% Copyright (C) 2011-2014  Jason McEwen
% See LICENSE.txt for license details

clear all;
close all;

% Define parameters.
L = 64;
el = 4
m = 2
gamma = pi/2
beta = pi/4
alpha = -pi/2
type = 'colour';

% Generate spherical harmonics.
flm = zeros(L^2,1);
ind = ssht_elm2ind(el, m);
flm(ind) = 1.0;

% Compute function on the sphere.
f = ssht_inverse(complex(real(flm), imag(flm)), L);

% Compute sampling grids.
%[thetas, phis, n, ntheta, nphi] = ssht_sampling(L, 'Grid', true);

% Precompute Wigner small-d functions
d = zeros(L, 2*L-1, 2*L-1);
d(1,:,:) = ssht_dl(squeeze(d(1,:,:)), L, 0, beta);
for el = 1:L-1
    d(el+1,:,:) = ssht_dl(squeeze(d(el,:,:)), L, el, beta);
end

% Rotate spherical harmonic
flm_rot = ssht_rotate_flm_dl(flm, d, alpha, gamma);

% Compute rotated function on the sphere.
f_rot = ssht_inverse(complex(real(flm_rot), imag(flm_rot)), L);

% Plot function on sphere.
figure(1);
subplot(2,2,1)
ssht_plot_sphere(real(f), L, 'Type', type, ...
    'ColourBar', false, 'Lighting', true);
subplot(2,2,2)
ssht_plot_sphere(imag(f), L, 'Type', type, ...
    'ColourBar', false, 'Lighting', true);
% Plot rotated function on sphere.
subplot(2,2,3)
ssht_plot_sphere(real(f_rot), L, 'Type', type, ...
    'ColourBar', false, 'Lighting', true);
subplot(2,2,4)
ssht_plot_sphere(imag(f_rot), L, 'Type', type, ...
    'ColourBar', false, 'Lighting', true);
