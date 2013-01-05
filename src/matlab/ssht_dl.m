function [dl] = ssht_dl(dl, L, el, theta, varargin)
% ssht_dl - Compute Wigner small-d functions
%
% Compute the el-th plane of the Wigner small-d functions
% (from the (el-1)-th plane) using Risbo's method.
%
% Default usage is given by
%
%   dl = ssht_dl(dl, L, el, theta, <options>)
%
% where dl is the Wigner plane for all m and n of size
% (2*L-1)*(2*L-1), L is the harmonic band-limit, el is the current
% harmonic degree (i.e. dl input should already be computed for el-1,
% and dl output will be computed for el), and theta is the argument of
% the Wigner small-d functions.
%
% Options consist of parameter type and value pairs.  Valid options
% include:
%  'SqrtTable'      = { Precomputed square-roots from 0 to 2*(L-1)+1 }
%  'SignTable'      = { Precomputed (-1)^m signs from m=0 to L-1 }
%
% Author: Jason McEwen (www.jasonmcewen.org)

% SSHT package to perform spin spherical harmonic transforms
% Copyright (C) 2013  Jason McEwen
% See LICENSE.txt for license details

% Parse arguments.
p = inputParser;
p.addRequired('dl', @isnumeric);
p.addRequired('L', @isnumeric);
p.addRequired('el', @isnumeric);
p.addRequired('theta', @isnumeric);
p.addParamValue('SqrtTable', 0, @isnumeric);
p.addParamValue('SignTable', 0, @isnumeric);
p.parse(dl, L, el, theta, varargin{:});
args = p.Results;

% Compute sqrt_tbl if not passed as argument.

sqrt_tbl = sqrt([0:2*(L-1)+1]);
signs = ones(L-1,1);
signs[2:2:end] = -1;

% Compute signs if not passed as argument.


% Compute Wigner small-d functions.
dl = ssht_dl_mex(dl, L, el, theta, sqrt_tbl, signs);
  
