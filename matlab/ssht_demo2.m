%ssht_demo2

clear all;
%close all;

L = 5
spin = 0
method = 'DH'
close_plot = true
plot_samples = false

[thetas, phis, n, ntheta, nphi] = ssht_sampling(L, 'Method', method, ...
                                                'Grid', true);

f = sin(thetas);

figure;
ssht_plot_sphere(f, L, 'Method', method, 'Close', close_plot, ...
                 'PlotSamples', plot_samples);
