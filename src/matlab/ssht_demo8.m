% ssht_demo8 - Run demo8
%
% Evaluate Wigner and spherical harmonic functions.
% 
% Default usage is given by
%
%   ssht_demo8
%
% Author: Jason McEwen (www.jasonmcewen.org)

% SSHT package to perform spin spherical harmonic transforms
% Copyright (C) 2013  Jason McEwen
% See LICENSE.txt for license details

clear;

% Define parameters.
L = 16
theta = rand * pi;

% Precompute data.
sqrt_tbl = sqrt([0:2*(L-1)+1])';
signs = ones(L+1,1);
signs(2:2:end) = -1; 

% Evaluate Wigner functions.
dl = zeros(2*L-1,2*L-1);
dl2 = zeros(L, 2*L-1);
dl2m1 = zeros(L, 2*L-1);
err_dl = zeros(L,1);
for el = 0:L-1
  
  % Evaluate Wigner functions using Risbo's method.
  dl = ssht_dl(dl, L, el, theta, ...
    'SqrtTable', sqrt_tbl, 'SignTable', signs);  
  
  % Evaluate Wigner functions using 3-term recursion (e.g. Kostelec).
  for mm = -el:el
    tmp = dl2(:,mm+L);
    dl2(:,mm+L) = ssht_dln(dl2(:,mm+L), dl2m1(:,mm+L), ...
      L, el, mm, theta, ...
      'SqrtTable', sqrt_tbl, 'SignTable', signs);
    dl2m1(:,mm+L) = tmp;
  end
  
  % Compute error between different recursions.
  err_dl(el+1) = max(max(abs(dl2(:,:) - dl(L:2*L-1,:))));
  
end

max_err_dl = max(err_dl)

% Compute sampling grids.
[thetas, phis, n, ntheta, nphi] = ssht_sampling(L, 'Grid', true);

% Compute all spherical harmonic using different recurisions.
ylm_matlab = ssht_ylm(L, thetas, phis, 'Recursion', 'Matlab');
ylm_risbo = ssht_ylm(L, thetas, phis, 'Recursion', 'Risbo');
ylm_kostelec = ssht_ylm(L, thetas, phis, 'Recursion', 'Kostelec');
ylm_nr = ssht_ylm(L, thetas, phis, 'Recursion', 'NumericalRecipes');

% Compare with spherical harmonics computed by inverse transform.
ylm_check = cell(L^2, 1);
err_ylm_matlab = zeros(L^2,1);
err_ylm_risbo = zeros(L^2,1);
err_ylm_kostelec = zeros(L^2,1);
err_ylm_nr = zeros(L^2,1);
for el = 0:L-1
  for m = 0:el

    % Generate spherical harmonics.
    flm = zeros(L^2,1);
    ind = ssht_elm2ind(el, m);
    flm(ind) = 1.0;

    % Compute function on the sphere.
    ylm_check{ind} = ssht_inverse(complex(real(flm), imag(flm)), L);

    % Compute errors.
    err_ylm_matlab(ind) = sqrt(sum(sum(abs(ylm_matlab{ind} - ylm_check{ind}).^2)));
    err_ylm_risbo(ind) = sqrt(sum(sum(abs(ylm_risbo{ind} - ylm_check{ind}).^2)));
    err_ylm_kostelec(ind) = sqrt(sum(sum(abs(ylm_kostelec{ind} - ylm_check{ind}).^2)));
    err_ylm_nr(ind) = sqrt(sum(sum(abs(ylm_nr{ind} - ylm_check{ind}).^2)));
    
  end
end

max_err_ylm_matlab = max(err_ylm_matlab)
max_err_ylm_risbo = max(err_ylm_risbo)
max_err_ylm_kostelec = max(err_ylm_kostelec)
max_err_ylm_nr = max(err_ylm_nr)
