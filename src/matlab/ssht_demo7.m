% ssht_demo7 - Run demo7
%
% Integrate a band-limited function on the sphere using the symmetrised
% quadrature weights.
% 
% Default usage is given by
%
%   ssht_demo7
%
% Author: Jason McEwen (www.jasonmcewen.org)

% SSHT package to perform spin spherical harmonic transforms
% Copyright (C) 2011  Jason McEwen
% See LICENSE.txt for license details

function ssht_demo7

clear all;

% Define parameters.
method = 'MW'
L = 64
reality = false
renorm_plot = true

% Plotting paramters.
line_width = 1.1;
line_width_thick = 2.5;
marker_size = 5;
marker_type1 = '^';
marker_type2 = 'd';
green_light = [0.2 0.6 0.4];% x339966
green_dark = [0 0.4 0.2];   % x006633
blue_light = [0.2 0.4 0.8]; % x3366CC
blue_dark = [0 0 1];        % x0000FF
red_light = [1 0.4 0.2];    % xFF6633
red_dark = [0.8 0.2 0];     % xCC3300
grey_light = [0.4 0.4 0.4]; % x666666
grey_dark = [0.2 0.2 0.2];  % x333333
magenta_light = [0.8 0.4 0.6]; % xCC6699
magenta_dark = [0.6 0.2 0.4];  % x993366

% Generate random flms (of complex signal).
flm = zeros(L^2,1);
flm = rand(size(flm)) + sqrt(-1)*rand(size(flm));
flm = 2.*(flm - (1+sqrt(-1))./2);

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

% Compute weights.
for m = -(L-1):(L-1)
   w(m+L) = mw_weights(m);
end

% Apply phase shift to account for theta offset.
wr = w .* exp(-i*[-(L-1):(L-1)]*pi/(2.*L-1));

% Compute weights as function of theta.
wr = real(fft(ifftshift(wr))).*2.*pi./(2*L-1).^2;
t=0:2*L-2;
theta = pi*(2*t+1)./(2.*L-1);
figure;
if renorm_plot
   plot(theta, [sin(theta(1:L)), zeros(1,L-1)], ...
      'Color', grey_dark, ...
      'LineStyle', '-', ...
      'LineWidth', line_width);
   hold on;
   plot(theta, wr./(2*pi).^2.*(2*L-1).^2, ...
      'Color', red_dark, ...
      'Marker', marker_type1, ...
      'MarkerSize', marker_size, ...
      'MarkerFaceColor', red_light, ...
      'MarkerEdgeColor', red_dark, ...
      'LineStyle', 'none', ...
      'LineWidth', line_width);
else
   plot(theta, (2*pi).^2/(2*L-1).^2.*[sin(theta(1:L)), zeros(1,L-1)],...
      'Color', grey_dark, ...      
      'LineStyle', '-', ...
      'LineWidth', line_width);      
   hold on;
   plot(theta, wr, ...      
      'Color', red_dark, ...
      'Marker', marker_type1, ...
      'MarkerSize', marker_size, ...
      'MarkerFaceColor', red_light, ...
      'MarkerEdgeColor', red_dark, ...
      'LineStyle', 'none', ...
      'LineWidth', line_width);

end

% Compute symmetrised quadrature weights defined on sphere.
q = wr(1:L);
q(1:L-1) = q(1:L-1) + wr(end:-1:L+1);
if renorm_plot
   plot(theta(1:L), q./(2*pi).^2.*(2*L-1).^2, ...
      'Color', magenta_dark, ...
      'Marker', marker_type2, ...
      'MarkerSize', marker_size, ...
      'MarkerFaceColor', magenta_light, ...
      'MarkerEdgeColor', magenta_dark, ...
      'LineStyle', 'none',  ...
      'LineWidth', line_width);
   a = axis; a(1) = 0; a(2) = 2*pi;
   axis(a)
   xlabel('t')
   ylabel('w')
   figure
   plot(theta(1:L), q./(2*pi).^2.*(2*L-1).^2-sin(theta(1:L)), ...
      'Color', grey_dark, ...
      'Marker', marker_type2, ...
      'MarkerSize', marker_size, ...
      'MarkerFaceColor', grey_light, ...
      'MarkerEdgeColor', grey_dark, ...
      'LineStyle', '-',  ...
      'LineWidth', line_width);   
   a = axis; a(1) = 0; a(2) = pi; a(3) = -a(4);
   axis(a)
   xlabel('t')
   ylabel('w')
else
   plot(theta(1:L), q, ...
      'Color', magenta_dark, ...
      'Marker', marker_type2, ...
      'MarkerSize', marker_size, ...
      'MarkerFaceColor', magenta_light, ...
      'MarkerEdgeColor', magenta_dark, ...
      'LineStyle', 'none',  ...
      'LineWidth', line_width);   
   a = axis; a(1) = 0; a(2) = 2*pi;
   axis(a)
   xlabel('t')
   ylabel('w')
   figure
   plot(theta(1:L), q - (2*pi).^2/(2*L-1).^2.*sin(theta(1:L)), ...
      'Color', grey_dark, ...
      'Marker', marker_type2, ...
      'MarkerSize', marker_size, ...
      'MarkerFaceColor', grey_light, ...
      'MarkerEdgeColor', grey_dark, ...
      'LineStyle', '-',  ...
      'LineWidth', line_width);    
   a = axis; a(1) = 0; a(2) = pi; a(3) = -a(4);
   axis(a)
   xlabel('t')
   ylabel('w')
end

% Integral of function given by rescaled (el,m)=(0,0) harmonic coefficient.
I0 = flm(1,1) .* sqrt(4*pi)

% Integrate function on sphere using all points.
f = ssht_inverse(flm, L, 'Method', method, 'Reality', reality);
Q = q.' * ones(1, 2*L-1);
I1 = sum(Q(:) .* f(:))

% Integration function on sphere using L points in theta.
% (Zero-pad to evalue the band-limited function on the correct L points.)
for r = 1:size(f,1)
   f_up(r,:) = ifft(ifftshift([0, fftshift(fft(f(r,:)))])) ...
      ./ (2*L-1) * (2*L);
end
f_down = f_up(:,1:2:end);
Q_down = q.' * ones(1, L);
I2 = sum(Q_down(:) .* f_down(:)) * (2*L-1) / L

% Compute integration errors.
I1_err = max(abs(I1 - I0))
I2_err = max(abs(I2 - I0))



function w = mw_weights(m)
%% mw_weights - Compute MW weights
%
% Compute the MW weights.
% 
% Default usage is given by
%
%   w = mw_weights(m)
%
% where m is the weight index from -(L-1):(L-1) and w is the corresponding
% weight.
%
% Author: Jason McEwen (www.jasonmcewen.org)

if m == 1
   w = sqrt(-1) * pi / 2;
elseif m == -1
   w = - sqrt(-1) * pi / 2;
elseif mod(m,2) == 0
   w = 2 / (1 - m.*m);
else
   w = 0;
end



