function ssht_plot_sphere(f, L, varargin)
% ssht_plot_sphere - Plot function on sphere
%
% Plots a functions defined on the sphere for various sampling schemes.
%
% Default usage is given by
%
%   ssht_plot_sphere(f, L, <options>)
%
% where f is the sampled function and L is the harmonic band-limit.
%
% Options consist of parameter type and value pairs.  Valid options
% include:
%  'Method'          = { 'MW'         [McEwen & Wiaux sampling (default)],
%                        'MWSS'       [McEwen & Wiaux symmetric sampling],
%                        'DH'         [Driscoll & Healy sampling],
%                        'GL'         [Gauss-Legendre sampling] }
%  'Type'            = { 'colour'     [generate colour plot (default)],
%                        'parametric' [generate parametric plot] }
%  'Close'           = { true         [close plot in phi (default)],
%                        false        [do not close plot in phi] }
%  'PlotSamples'     = { false        [do not plot sample positions (default)],
%                        true         [plot sample positions] }
%  'ParametricScale' = scale          [scaling for parametric plot (default=0.5)]
%  'Lighting'        = { false        [do not light plot (default)],
%                        true         [light plot] }

% Parse arguments.
p = inputParser;
p.addRequired('f', @isnumeric);     
p.addRequired('L', @isnumeric);          
p.addParamValue('Type', 'colour', ...
   @(x)strcmpi(x,'colour') | strcmpi(x,'parametric'));
p.addParamValue('Method', 'MW', @ischar);
p.addParamValue('Close', true, @islogical);
p.addParamValue('PlotSamples', false, @islogical);
p.addParamValue('ParametricScale', 0.5, @isnumeric);
p.addParamValue('Lighting', false, @islogical);
p.parse(f, L, varargin{:});
args = p.Results;

% Define parameters.
TOL = 1e-10;
PARAMETRIC_SCALE = args.ParametricScale;

% Compute grids.
minf = min(f(:));
maxf = max(f(:));
[thetas, phis, n, ntheta, nphi] = ssht_sampling(L, 'Method', ...
                                                args.Method, 'Grid', true);
if (size(thetas) ~= size(f)) 
  error('Size of f does not match sampling of method %s.', args.Method);
end

% Compute position scaling for parametric plot.
if strcmpi(args.Type, 'parametric')
   if abs(maxf - minf) < TOL
      f_normalised = f;
   else
      f_normalised = (f - minf)./(maxf - minf).*PARAMETRIC_SCALE + 1.0;
   end
else
   f_normalised = ones(size(f));
end

% Close plot.
if args.Close
   close = @(x) [x, x(:,1)];
   f = close(f);
   f_normalised = close(f_normalised);
   thetas = close(thetas);
   phis = close(phis);
end

% Compute location of vertices.
[x, y, z] = ssht_coord_s2c(thetas, phis, 1.0);
x = x .* f_normalised;
y = y .* f_normalised;
z = z .* f_normalised;

% Plot.
h = surf(x,y,z,f);
caxis([minf, maxf]);
%colorbar('vert');
set(h, 'LineStyle', 'none')
if args.PlotSamples
    hold on;
    plot3(x, y, z, 'k.');
    hold off;
end
axis([-1,1,-1,1,-1,1]);
view([-1,1,1]);
axis equal;
if args.Lighting
  axis off; 
  light; 
  lighting phong; 
  camzoom(1.3);
end


