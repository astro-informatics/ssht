function [flm] = ssht_forward(f, L, varargin)
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
p.addRequired('f', @isnumeric);          
p.addRequired('L', @isnumeric);          
p.addParamValue('Method', 'MW', @ischar);
p.parse(f, L, varargin{:});
args = p.Results;

% Computing sampling points.

spin = 0
verbosity = 0

%f = ssht_inverse_mex(flm, L, args.Method, spin, verbosity);

flm = ssht_forward_mex(f, L, args.Method, spin, verbosity);




  