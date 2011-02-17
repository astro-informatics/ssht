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
%
% Author: Jason McEwen (www.jasonmcewen.org)

% Parse arguments.
p = inputParser;
p.addRequired('f', @isnumeric);     
p.addRequired('L', @isnumeric);          
p.addParamValue('Method', 'MW', @ischar);
p.addParamValue('ColourBar', false, @islogical);
p.parse(f, L, varargin{:});
args = p.Results;

% Compute grids.
minf = min(f(:));
maxf = max(f(:));
[thetas, phis, n, ntheta, nphi] = ssht_sampling(L, 'Method', ...
                                                args.Method, 'Grid', true);
if (size(thetas) ~= size(f)) 
  error('Size of f does not match sampling of method %s.', args.Method);
end

% Compute Mollweide projection.
[x, y] = ssht_mollweide(thetas, phis);

% Plot.
h = surf(x,y,zeros(size(x)),f);
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
