%ssht_demo2

clear all;
%close all;

L = 5
spin = 0
method = 'MW'
close_plot = true
plot_samples = false
reality = true

[thetas, phis, n, ntheta, nphi] = ssht_sampling(L, 'Method', method, ...
                                                'Grid', true);

f = sin(thetas);

figure;
ssht_plot_sphere(f, L, 'Method', method, 'Close', close_plot, ...
                 'PlotSamples', plot_samples, 'Lighting', true);


% Random flms (of complex signal).
flm = zeros(L^2,1);
flm = rand(size(flm)) + sqrt(-1)*rand(size(flm));
flm = 2.*(flm - (1+sqrt(-1))./2);
ind_min = spin^2 + abs(spin);
flm(1:ind_min) = 0;

% Impose reality on flms.
if reality
   for el = 0:L-1
      m = 0;
      ind = el^2 + el + m;
      flm(ind+1) = real(flm(ind+1));
      for m = 1:el
         ind_pm = el^2 + el + m;
         ind_nm = el^2 + el - m;
         flm(ind_nm+1) = (-1)^m * conj(flm(ind_pm+1));
      end
   end
end

f = ssht_inverse(flm, L, 'Method', method);