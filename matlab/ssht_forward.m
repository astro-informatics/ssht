function [flm] = ssht_forward(f, L, varargin)
% ssht_forward - Compute forward spin spherical harmonic transform
%
% Computes forward spin spherical harmonic transform for various
% exact sampling theorems.
%
% Default usage is given by
%
%   flm = ssht_forward(f, L, <options>)
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
p.addRequired('f', @isnumeric);          
p.addRequired('L', @isnumeric);          
p.addParamValue('Method', 'MW', @ischar);
p.addParamValue('Spin', 0, @isnumeric);
p.addParamValue('Reality', false, @islogical);
p.parse(f, L, varargin{:});
args = p.Results;

% Computing forward transform.
flm = ssht_forward_mex(f, L, args.Method, args.Spin, args.Reality);




  