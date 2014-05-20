function ssht_plot_mollweide(f, L, varargin)
% ssht_plot_mollweide - Plot function using Mollweide projection
%
% Plots a Mollweide projection of a function defined on the sphere, for
% various sampling schemes. 
%
% Default usage is given by
%
%   ssht_plot_mollweide(f, L, <options>)
%
% where f is the sampled function and L is the harmonic band-limit.
%
% Options consist of parameter type and value pairs.  Valid options
% include:
%  'Method'          = { 'MW'         [McEwen & Wiaux sampling (default)],
%                        'MWSS'       [McEwen & Wiaux symmetric sampling],
%                        'DH'         [Driscoll & Healy sampling],
%                        'GL'         [Gauss-Legendre sampling] }
%  'ColourBar'       = { false        [do not add colour bar (default)],
%                        true         [add colour bar] }
%  'Mode'            = { 0            Plot amplitude, or modulus is f complex (default),
%                        1            Plot real part,
%                        2            Plot imaginaty part,
%                        3            Plot modulus and arrows for real/img angle }
%  'Spin'            = { non-negative integers (default=0) }
%
% Author: Jason McEwen (www.jasonmcewen.org)

% SSHT package to perform spin spherical harmonic transforms
% Copyright (C) 2011  Jason McEwen
% See LICENSE.txt for license details

% Parse arguments.
p = inputParser;
p.addRequired('f', @isnumeric);     
p.addRequired('L', @isnumeric);          
p.addParamValue('Method', 'MW', @ischar);
p.addParamValue('ColourBar', false, @islogical);
p.addParamValue('Mode', 0, @isnumeric);
p.addParamValue('SubL', 16, @isnumeric);
p.parse(f, L, varargin{:});
args = p.Results;

SubL = min([args.SubL L]);

% Compute grids.
minf = min(f(:));
maxf = max(f(:));
[thetas, phis, n, ntheta, nphi] = ssht_sampling(L, 'Method', args.Method, 'Grid', true);
if (size(thetas) ~= size(f)) 
  error('Size of f does not match sampling of method %s.', args.Method);
end

% Compute Mollweide projection.
[x, y] = ssht_mollweide(thetas, phis);

switch args.Mode
    case {1}
        op = str2func('real');
    case {2}
        op = str2func('imag');
    otherwise
        op = str2func('abs');
end

% Plot.
h = surf(x,y,zeros(size(x)),op(f));

if args.Mode == 3
    hold on
    f_lm = ssht_forward(f, L, 'spin', 2, 'Method', args.Method);
    subf_lm = f_lm(1:SubL^2.0,1);
    subf = ssht_inverse(subf_lm, SubL, 'spin', 2, 'Method', 'GL');
    [thetas, phis, n, ntheta, nphi] = ssht_sampling(SubL, 'Method', 'GL', 'Grid', true); 
    [subx, suby] = ssht_mollweide(thetas, phis);  
    quiver(subx,suby,real(subf),imag(subf),'black', 'AutoScaleFactor', 1.0, 'ShowArrowHead', 'off', 'LineWidth', 1);
end    

set(h, 'LineStyle', 'none')
axis equal
axis off
campos([0 0 1])
camup([0 1 0])


function [x, y] = ssht_mollweide(thetas, phis)
% ssht_mollweide - Compute Mollweide projection
%
% Compute Mollweide projection of spherical coordinates.
%
% Usage is given by
%
%   [x,y] = ssht_mollweide(thetas, phis)
%
% where thetas and phis are spherical coordinates and x and y are the
% projected Mollweide coordinates.
%
% Author: Jason McEwen (www.jasonmcewen.org)

MAX_ITERATIONS = 1e5;
TOL = 1e-10;

% Convert theta to longitude.
thetas = pi/2 - thetas;
phis = phis - pi;

t = thetas;
for it = 1:MAX_ITERATIONS

   dt = (t + sin(t) - pi.*sin(thetas)) ./ (1 + cos(t));
   t = t - dt;
   
   if(max(abs(dt)) < TOL)
      break;
   end
   
end
t = t/2;
x = 2 .* sqrt(2) ./ pi .* phis .* cos(t);
y = sqrt(2) .* sin(t);
