function [thetas, phis, varargout] = ssht_sampling(L, varargin)
% ssht_sampling - Compute sample positions
%
% Computes sample positions on sphere for various sampling methods.
%
% Default usage is given by
%
%   [thetas, phis] = ssht_sampling(L, <options>)
%
% where L is the harmonic band-limit and theta and phi specify sample
% positions. 
%
% Options consist of parameter type and value pairs.  Valid options
% include:
%  'Method' = { 'MW'   [McEwen & Wiaux sampling (default)],
%               'MWSS' [McEwen & Wiaux symmetric sampling],
%               'DH'   [Driscoll & Healy sampling],
%               'GL'   [Gauss-Legendre sampling] }
%  'Grid'   = { false  [return theta and phi vectors (default)],
%               true   [return theta and phi grids] }
% 
% May optionally return total number of samples n, number of theta
% samples ntheta and number of phi samples nphi through usage
%
%   [thetas, phis, n, ntheta, nphi] = ssht_sampling(L, <options>)
%
% Author: Jason McEwen (jason.mcewen@epfl.ch)

% Parse arguments.
p = inputParser;
p.addRequired('L', @isnumeric);          
p.addParamValue('Method', 'MW', @ischar);
p.addParamValue('Grid', false, @islogical);
p.parse(L, varargin{:});
args = p.Results;

% Computing sampling points.
[n, ntheta, nphi, thetas, phis] = ssht_sampling_mex(L, args.Method);

% Compute grids if requested.
if args.Grid
  [thetas, phis] = ndgrid(thetas, phis);
end

% Set optional outputs.
if nargout >= 3
  varargout(1) = {n};
end
if nargout >= 4
  varargout(2) = {ntheta};
end
if nargout >= 5
  varargout(3) = {nphi};
end

  