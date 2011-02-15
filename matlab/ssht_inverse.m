function [f] = ssht_inverse(flm, L, varargin)
% ssht_inverse - Compute inverse spin spherical harmonic transform
%
% Computes inverse spin spherical harmonic transform for various
% exact sampling theorems.
%
% Default usage is given by
%
%   f = ssht_inverse(flm, L, <options>)
%
% where L is the harmonic band-limit, flm is the vector of L^2
% harmonic coefficients and f is the sampled function values
% indexed by theta and phi.
%
% Options consist of parameter type and value pairs.  Valid options
% include:
%  'Method'          = { 'MW'         [McEwen & Wiaux sampling (default)],
%                        'MWSS'       [McEwen & Wiaux symmetric sampling],
%                        'DH'         [Driscoll & Healy sampling],
%                        'GL'         [Gauss-Legendre sampling] }
%  'Spin'            = { non-negative integers (default=0) }
%  'Reality'         = { false        [do not assume f real (default)],
%                        true         [assume f real (improves performance)] }
% 
% Author: Jason McEwen (jason.mcewen@epfl.ch)

% Parse arguments.
p = inputParser;
p.addRequired('flm', @isnumeric);          
p.addRequired('L', @isnumeric);          
p.addParamValue('Method', 'MW', @ischar);
p.addParamValue('Spin', 0, @isnumeric);
p.addParamValue('Reality', false, @islogical);
p.parse(flm, L, varargin{:});
args = p.Results

% Computing inverse transform.
f = ssht_inverse_mex(flm, L, args.Method, args.Spin, args.Reality);




  