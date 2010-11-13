function ssht_demo
%%
% SSHT_DEMO - Compute forward and inverse spherical harmonic transforms
% 
% Harmonic transforms are computed by calling the fortran programs
% ssht_forward and ssht_inverse through the command line (these 
% executables should be contained in the same directory as ssht_demo.m).
%
% Author: Jason McEwen (jason.mcewen@epfl.ch)
% Date: November 2010


% Define size parameters.
spin = -4;
L = 16;
nphi = 2*L - 1;
ntheta_dh = 2*L;
ntheta_mw = L;

% Define sample points.
t_dh = 0:ntheta_dh-1;
theta_dh =  (2d0*t_dh+1d0)*pi / (4d0*L)
t_mw = 0:ntheta_mw-1;
theta_mw = (2d0*t_mw+1d0)*pi / (2d0*L - 1d0)
p = 0:nphi-1;
phi_dh = 2d0*p*pi / (2d0*L - 1d0)
phi_mw = (2d0*p+1d0)*pi / (2d0*L - 1d0)

% Random flms.
flm = zeros(L^2,1);
flm = rand(size(flm)) + sqrt(-1)*rand(size(flm));
flm = 2.*(flm - (1+sqrt(-1))./2);
ind_min = spin^2 + abs(spin);
flm(1:ind_min) = 0;

% Measure error of inverse-forward transform for DH.
f_dh = ssht_inverse(flm, 'DH', L, spin);
flm_dh_syn = ssht_forward(f_dh, 'DH', L, spin);
maxerr_dh = max(abs(flm - flm_dh_syn))

% Measure error of inverse-forward transform for MW.
f_mw = ssht_inverse(flm, 'MW', L, spin);
f_mw_grid = reshape(f_mw, nphi, ntheta_mw);
f_mw_grid(:,end) = f_mw_grid(1,end).*exp(sqrt(-1)*spin*phi_mw).* ...
   exp(sqrt(-1)*(-1d0)*pi / (2d0*L - 1d0)*spin);
f_mw = f_mw_grid(:);
flm_mw_syn = ssht_forward(f_mw, 'MW', L, spin);
maxerr_mw = max(abs(flm - flm_mw_syn))

% Define sample points and function values on grids.
f_dh_grid = reshape(f_dh, nphi, ntheta_dh)
[theta_dh_grid, phi_dh_grid] = meshgrid(theta_dh, phi_dh)
f_mw_grid = reshape(f_mw, nphi, ntheta_mw)
[theta_mw_grid, phi_mw_grid] = meshgrid(theta_mw, phi_mw)

keyboard;


function flm = ssht_forward(f, method, L, spin)
%% Compute forward spherical harmonic transform for method={'DH','MW'}.

f_save = [real(f), imag(f)];
save('f.txt', 'f_save', '-ascii', '-double');

cmd = sprintf('%s%s%s%s%s%s', ...
   '../bin/ssht_forward -inp f.txt -out flm.txt -method ', ...
   method(1:2), ' -L ', num2str(L), ' -spin ', num2str(spin));
system(cmd);

flm = load('flm.txt');
flm = flm(:,1) + sqrt(-1)*flm(:,2);

system('rm -f flm.txt f.txt');


function f = ssht_inverse(flm, method, L, spin)
%% Compute inverse spherical harmonic transform for method={'DH','MW'}.

flm_save = [real(flm), imag(flm)];
save('flm.txt', 'flm_save', '-ascii', '-double');

cmd = sprintf('%s%s%s%s%s%s', ...
   '../bin/ssht_inverse -inp flm.txt -out f.txt -method ', ...
   method(1:2), ' -L ', num2str(L), ' -spin ', num2str(spin));
system(cmd);

f = load('f.txt');
f = f(:,1) + sqrt(-1)*f(:,2);

system('rm -f flm.txt f.txt');