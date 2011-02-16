% ssht_demo4 - Run demo4
%
% Demo to compute inverse and forward transform of spin function, using
% polar interface with various options.  
%
% Default usage is given by
%
%   ssht_demo4
%
% Author: Jason McEwen (jason.mcewen@epfl.ch)

clear all;

% Define parameters.
L = 64
spin = 0
reality = true
methods = {'MW', 'MWSS'};
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
if strcmp(method, 'MW')
   [f, f_sp, phi_sp] = ssht_inverse(flm, L, ...
      'Method', method, ...
      'Spin', spin, ...
      'Reality', reality);
   flm_syn = ssht_forward(f, L, ...
      'Method', method, ...
      'Spin', spin, ...
      'Reality', reality, ...
      'SouthPoleSample', f_sp, ...
      'SouthPolePhi', phi_sp);
else
   [f, f_sp, phi_sp, f_np, phi_np] = ssht_inverse(flm, L, ...
      'Method', method, ...
      'Spin', spin, ...
      'Reality', reality);
   flm_syn = ssht_forward(f, L, ...
      'Method', method, ...
      'Spin', spin, ...
      'Reality', reality, ...
      'SouthPoleSample', f_sp, ...
      'SouthPolePhi', phi_sp, ...
      'NorthPoleSample', f_np, ...
      'NorthPolePhi', phi_np);
end         
                            
% Compute max error in harmonic space.
maxerr = max(abs(flm_syn - flm))

