function [dlp1] = ssht_dln(dl, dlm1, L, el, n, theta, varargin)
% ssht_dln - Compute Wigner small-d functions for given n
%
% Compute the el-th line of the Wigner small-d functions for given n
% (from the (el-1)-th and (el-2)-th lines) using 3-term recursion of 
% Kostelec.
%
% Default usage is given by
%
%   dlp1 = ssht_dln(dl, dlm1, L, el, n, theta, <options>)
%
% where dl is the Wigner line for el for non-negative m and given n of size
% L, dlm1 is the line for el-1 and dlp1 is the line computed for el+1,L is
% the harmonic band-limit, el is the current harmonic degree, and theta is 
% the argument of the Wigner small-d functions.
%
% Options consist of parameter type and value pairs.  Valid options
% include:
%  'SqrtTable'      = { Precomputed square-roots from 0 to 2*(L-1)+1 }
%  'SignTable'      = { Precomputed (-1)^m signs from m=0 to L }
%
% Author: Jason McEwen (www.jasonmcewen.org)

% SSHT package to perform spin spherical harmonic transforms
% Copyright (C) 2013  Jason McEwen
% See LICENSE.txt for license details

% Parse arguments.
p = inputParser;
p.addRequired('dl', @isnumeric);
p.addRequired('dlm1', @isnumeric);
p.addRequired('L', @isnumeric);
p.addRequired('el', @isnumeric);
p.addRequired('n', @isnumeric);
p.addRequired('theta', @isnumeric);
p.addParamValue('SqrtTable', 0, @isnumeric);
p.addParamValue('SignTable', 0, @isnumeric);
p.parse(dl, dlm1, L, el, n, theta, varargin{:});
args = p.Results;
defaults = p.UsingDefaults;

% Compute sqrt_tbl if not passed as argument.
if (sum(ismember(defaults, 'SqrtTable')))
  sqrt_tbl = sqrt([0:2*(L-1)+1])';
else
  sqrt_tbl = args.SqrtTable;
end

% Compute signs if not passed as argument.
if (sum(ismember(defaults, 'SqrtTable')))
  signs = ones(L+1,1);
  signs(2:2:end) = -1; 
else
  signs = args.SignTable;
end

% Compute Wigner small-d functions.
dlp1 = ssht_dln_mex(dl, dlm1, L, el, n, theta, sqrt_tbl, signs);
