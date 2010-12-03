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
L = 12;
spin = 2;
reality = 0;
verbosity = 1;

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
phi_mw = 2d0*p*pi / (2d0*L - 1d0)

% Random flms (of complex signal).
flm = zeros(L^2,1);
flm = rand(size(flm)) + sqrt(-1)*rand(size(flm));
flm = 2.*(flm - (1+sqrt(-1))./2);
ind_min = spin^2 + abs(spin);
flm(1:ind_min) = 0;

% Impose reality on flms.
if (reality == 1)
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

% Measure error of inverse-forward transform for DH.
f_dh = ssht_inverse(flm, 'DH', L, spin, reality, verbosity);
f_sp = 0; phi_sp = 0;
flm_dh_syn = ssht_forward(f_dh, f_sp, phi_sp, 'DH', L, spin, reality, verbosity);
maxerr_dh = max(abs(flm - flm_dh_syn))

% Measure error of inverse-forward transform for GL.
% TODO...



% Measure error of inverse-forward transform for MW.
[f_mw, f_sp, phi_sp] = ssht_inverse(flm, 'MW', L, spin, reality, verbosity);
flm_mw_syn = ssht_forward(f_mw, f_sp, phi_sp, 'MW', L, spin, reality, verbosity);
maxerr_mw = max(abs(flm - flm_mw_syn))

% Define sample points and function values on grids.
f_dh_grid = reshape(f_dh, nphi, ntheta_dh);
[theta_dh_grid, phi_dh_grid] = meshgrid(theta_dh, phi_dh);
f_mw_grid = reshape(f_mw, nphi, ntheta_mw-1);
[theta_mw_grid, phi_mw_grid] = meshgrid(theta_mw(1:end), phi_mw);
%[theta_mw_grid, phi_mw_grid] = meshgrid(theta_mw(1:end-1), phi_mw);


% Define Gauss-Legendre sample points.
addpath('~/AP/src/chebfun-3.111')
phi_gl = 2d0*p*pi / (2d0*L - 1d0)
theta_gl = acos(legpts(L)); theta_gl = theta_gl(end:-1:1)';
[theta_gl_grid, phi_gl_grid] = meshgrid(theta_gl, phi_gl);


%% Plot sample points on sphere.

[x_dh, y_dh, z_dh] = ssht_s2c(theta_dh_grid, phi_dh_grid, 1.0);
[x_gl, y_gl, z_gl] = ssht_s2c(theta_gl_grid, phi_gl_grid, 1.0);
[x_mw, y_mw, z_mw] = ssht_s2c(theta_mw_grid, phi_mw_grid, 1.0);

figure;
r = 0.99;
surf(r.*x_dh, r.*y_dh, r.*z_dh, ones(size(x_dh)), 'EdgeColor', 'none');
hold on;

green_light = [0.2 0.6 0.4];% x339966
green_dark = [0 0.4 0.2];   % x006633
blue_light = [0.2 0.4 0.8]; % x3366CC
blue_dark = [0 0 1];        % x0000FF
red_light = [1 0.4 0.2];    % xFF6633
red_dark = [0.8 0.2 0];     % xCC3300

marker_size = 4;
plot3(x_dh, y_dh, z_dh, 'o', ...
   'MarkerSize', marker_size, ...
   'MarkerFaceColor', green_light, ...
   'MarkerEdgeColor', green_dark);
plot3(x_gl, y_gl, z_gl, 'o', ...
   'MarkerSize', marker_size, ...
   'MarkerFaceColor', blue_light, ...
   'MarkerEdgeColor', blue_dark);
plot3(x_mw, y_mw, z_mw, 'o', ...
   'MarkerSize', marker_size, ...
   'MarkerFaceColor', red_light, ...
   'MarkerEdgeColor', red_dark);

axis square;
axis off;
colormap('gray');
caxis([-1.0, 1.2]);

campos([-1 1 -0.8])
zoom(1.5)

% Plot longitude rings.
eq_phis = [0:0.1:2*pi];
for t = 0:pi/6:pi
   eq_thetas = ones(size(eq_phis))*t
   [eq_x, eq_y, eq_z] = ssht_s2c(eq_thetas, eq_phis, 1.0);
   plot3(eq_x, eq_y, eq_z, '-k')
end

% Plot lattitude rings.
eq_thetas = [0:0.1:pi];
for t = 0:pi/6:2*pi
   eq_phis = ones(size(eq_thetas))*t
   [eq_x, eq_y, eq_z] = ssht_s2c(eq_thetas, eq_phis, 1.0);
   plot3(eq_x, eq_y, eq_z, '-k')
end


%% Plot number of points required vs band-limit. 

line_width = 1.8;

fig1 = figure;
L = [1:1:8];
N_dh = 2.*L .* (2.*L - 1);
N_mw = (L-1) .* (2.*L - 1) + 1;
N_gl = L .* (2.*L - 1);
plot(L, N_dh, 'Color', green_dark, ...
   'LineWidth', line_width, 'LineStyle', '--'); 
hold on;
plot(L, N_gl, 'Color', blue_dark, ...
   'LineWidth', line_width, 'LineStyle', '-.');
plot(L, N_mw, 'Color', red_dark, ...
   'LineWidth', line_width, 'LineStyle', '-')
axis([L(1) L(end) 0 250])

fig2 = figure;
L = [1:1:4096];
N_dh = 2.*L .* (2.*L - 1);
N_mw = (L-1) .* (2.*L - 1) + 1;
N_gl = L .* (2.*L - 1);
plot(L, N_dh, 'Color', green_dark, ...
   'LineWidth', line_width, 'LineStyle', '--');
hold on;
plot(L, N_gl, 'Color', blue_dark, ...
   'LineWidth', line_width, 'LineStyle', '-.');
plot(L, N_mw, 'Color', red_dark, ...
   'LineWidth', line_width, 'LineStyle', '-')
axis([L(1) L(end) 0 7e7])
xlabel('L');
ylabel('N');

[h_m h_i]=inset(fig2, fig1, 0.5);


function [x, y, z] = ssht_s2c(theta, phi, r)
%% Compute cartesian corrdinates from spherical coordinates.

x = r .* sin(theta) .* cos(phi);
y = r .* sin(theta) .* sin(phi);
z = r .* cos(theta);


function flm = ssht_forward(f, f_sp, phi_sp, method, L, spin, reality, verbosity)
%% Compute forward spherical harmonic transform for method={'DH','MW'}.

d = [real(f), imag(f)];
d = [real(f_sp), imag(f_sp); d];
d = [real(phi_sp), imag(phi_sp); d];
save('f.txt', 'd', '-ascii', '-double');

cmd = sprintf('%s%s%s%s%s%s', ...
   '../bin/ssht_forward -inp f.txt -out flm.txt -method ', ...
   method(1:2), ' -L ', num2str(L), ' -spin ', num2str(spin), ...
   ' -reality ', num2str(reality), ' -verbosity ', num2str(verbosity));
system(cmd);

flm = load('flm.txt');
flm = flm(:,1) + sqrt(-1)*flm(:,2);

system('rm -f flm.txt f.txt');


function [f, f_sp, phi_sp] = ssht_inverse(flm, method, L, spin, reality, verbosity)
%% Compute inverse spherical harmonic transform for method={'DH','MW'}.

flm_save = [real(flm), imag(flm)];
save('flm.txt', 'flm_save', '-ascii', '-double');

cmd = sprintf('%s%s%s%s%s%s', ...
   '../bin/ssht_inverse -inp flm.txt -out f.txt -method ', ...
   method(1:2), ' -L ', num2str(L), ' -spin ', num2str(spin), ...
   ' -reality ', num2str(reality), ' -verbosity ', num2str(verbosity));
system(cmd);

d = load('f.txt');
phi_sp = d(1,1);
f = d(3:end,1);
f_sp = d(2,1);
if(reality ~= 1)
   f_sp = f_sp + sqrt(-1)*d(2,2);
   f = f + sqrt(-1)*d(3:end,2);
end

system('rm -f flm.txt f.txt');
