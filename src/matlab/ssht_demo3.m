% ssht_demo3 - Run demo3
%
% Demo to compute inverse and forward transform of spin function, using
% standard interface with various options.  
%
% Default usage is given by
%
%   ssht_demo3
%
% Author: Jason McEwen (www.jasonmcewen.org)

clear all;

% Define parameters.
L = 64
spin = 2
reality = false
methods = {'MW', 'MWSS', 'GL', 'DH'};
method = char(methods(1))
close_plot = true;
plot_samples = false;

% Generate random flms (of complex signal).
flm = zeros(L^2,1);
flm = rand(size(flm)) + sqrt(-1)*rand(size(flm));
flm = 2.*(flm - (1+sqrt(-1))./2);

% Zero harmonic coefficients with el<|spin|.
ind_min = spin^2 + abs(spin);
flm(1:ind_min) = 0;

% Impose reality on flms.
if reality
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
end

% Compute inverse then forward transform.
f = ssht_inverse(flm, L, 'Method', method, 'Spin', spin, ...
                 'Reality', reality);
flm_syn = ssht_forward(f, L, 'Method', method, 'Spin', spin, ...
                 'Reality', reality);

% Compute max error in harmonic space.
maxerr = max(abs(flm_syn - flm))

% Compute sampling grids.
[thetas, phis, n, ntheta, nphi] = ssht_sampling(L, 'Method', method, ...
                                                'Grid', true);

% Plot function on sphere.
figure;
ssht_plot_sphere(abs(f), L, 'Method', method, 'Close', close_plot, ...
                 'PlotSamples', plot_samples, 'Lighting', true);
