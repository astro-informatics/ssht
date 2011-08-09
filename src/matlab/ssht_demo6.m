% ssht_demo6 - Run demo6
%
% Smooth Earth topography map by applying a Gaussian filter in harmonic
% space.  
%
% The official Earth Gravitational Model EGM2008 has been publicly
% released by the U.S. National Geospatial-Intelligence Agency (NGA) (these
% data were downloaded from Frederik Simons' web page:
% http://geoweb.princeton.edu/people/simons).
% 
% Default usage is given by
%
%   ssht_demo6
%
% Author: Jason McEwen (www.jasonmcewen.org)

% SSHT package to perform spin spherical harmonic transforms
% Copyright (C) 2011  Jason McEwen
% See LICENSE.txt for license details

clear all;

% Define parameters.
methods = {'MW', 'MWSS', 'GL', 'DH'};
method = char(methods(1))
close_plot = true;
plot_samples = false;
sigma = 5e-4;

% Load harmonic coefficients of Earth.
load('data/EGM2008_Topography_flms_L0128');

% Smooth harmonic coefficients.
flm_smooth = flm;
for el = 1:L-1
   for m = 0:el
      ind = ssht_elm2ind(el,m);
      flm_smooth(ind) = flm_smooth(ind) ...
         .* exp(-el.^2 * sigma);
   end
end

% Compute real space version of Earth.
f = ssht_inverse(flm, L, 'Method', method, 'Reality', true);
f_smooth = ssht_inverse(flm_smooth, L, 'Method', method, 'Reality', true);
f = fliplr(f); f_smooth = fliplr(f_smooth);

% Plot functions on sphere.
figure;
ssht_plot_sphere(f, L, 'Method', method, 'Close', close_plot, ...
                 'PlotSamples', plot_samples, 'Lighting', true);
figure;
ssht_plot_mollweide(f, L, 'Method', method);
figure;
ssht_plot_sphere(f_smooth, L, 'Method', method, 'Close', close_plot, ...
                 'PlotSamples', plot_samples, 'Lighting', true);
figure;
ssht_plot_mollweide(f_smooth, L, 'Method', method);
