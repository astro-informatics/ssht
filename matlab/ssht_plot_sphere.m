% ******************************************************************************
% Function:     plot_sphere
% Description:  Plot the _magnitude_ of a function defined over the sphere.
% Notes:        f is the matrix of the function values evaluated over all 
%               combinations of theta and phi (obtained using meshgrid for 
%               example).
%               plot_type = {'parametric', 'colour'} selects appropriate type 
%               of plot.
%               plot_samples_posns = {'yes', 'no'}.
%               The argument plot_axes may not be included and new figures will
%               be generated.
% Author:       Jason McEwen
% Version:      1.01 (Last modified 13/02/03)
% ******************************************************************************


%
% 'Type' {'colour','parametric'}
% 'Title' string
% 'Close' true/false [default false]
% 'PlotSamples' true/false [default false]


function s2ea_plot_sphere(f, varargin)

% Parse arguments.
p = inputParser;
p.addRequired('f', @isnumeric);          
p.addParamValue('Type', 'colour', ...
   @(x)strcmpi(x,'colour') | strcmpi(x,'parametric'));
p.addParamValue('Close', false, @islogical);
p.addParamValue('PlotSamples', false, @islogical);
p.addParamValue('Title', '', @ischar);
p.addParamValue('ParametricScale', 0.5, @isnumeric);
p.parse(f, varargin{:});
args = p.Results;

% Define parameters.
TOL = 1e-10;
PARAMETRIC_SCALE = args.ParametricScale;

% Compute grids.
minf = min(f(:));
maxf = max(f(:));
B = s2ea_samples_f2b(f);
[theta_grid, phi_grid] = s2ea_samples_sgrid(B);

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
   close = @(x) [x; x(1,:)];
   f = close(f);
   f_normalised = close(f_normalised);
   theta_grid = close(theta_grid);
   phi_grid = close(phi_grid);
end

% Compute location of vertices.
[x, y, z] = s2ea_coord_s2c(theta_grid, phi_grid, 1.0);
x = x .* f_normalised;
y = y .* f_normalised;
z = z .* f_normalised;

% Plot.
h = surf(x,y,z,f);
caxis([minf, maxf]);
colorbar('vert');
set(h, 'LineStyle', 'none')
if args.PlotSamples
    hold on;
    plot3(x, y, z, 'k.');
    hold off;
end
axis([-1,1,-1,1,-1,1]);
view([-1,1,1]);
axis equal;
%axis off; 
%light; lighting phong; %camzoom(1.3);

if ~strcmp(args.Title, '')
   title(args.Title);
end

