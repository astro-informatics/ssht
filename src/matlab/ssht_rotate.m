function [f_rotated] = ssht_rotate(f, L, alpha, beta, gamma, varargin)
% ssht_rotate - Rotate a function defined on space
%
% Rotates a function on the sphere, where the rotation is specified by a 
% point (alpha, beta, gamma) in SO(3).
%
% Default usage is given by
%
%   [flm_rotated] = ssht_rotate(f, L, alpha, beta, gamma, <options>)
%
% where f is the function to be rotated, L is the harmonic band-limit, 
% alpha, beta and gamma are the Euler angles specifying the rotation, 
% and f_rotated the resulting function after rotation. 
%
% Authors: Jason McEwen (www.jasonmcewen.org)

% SSHT package to perform spin spherical harmonic transforms
% Copyright (C) 2011-2017 Jason McEwen
% See LICENSE.txt for license details

% Parse arguments.
p = inputParser;
p.addRequired('f', @isnumeric);          
p.addRequired('L', @isnumeric);    
p.addRequired('alpha', @isnumeric);          
p.addRequired('beta', @isnumeric);          
p.addRequired('gamma', @isnumeric);  
p.addParamValue('Method', 'MW', @ischar);
p.addParamValue('Spin', 0, @isnumeric);
p.addParamValue('Reality', false, @islogical);
p.parse(f, L, alpha, beta, gamma, varargin{:});
args = p.Results;

% Compute forward harmonic transform.
flm = ssht_forward(f, L, ...
  'Method', args.Method, 'Spin', args.Spin, 'Reality', args.Reality);

% Rotate in harmonic space.
flm_rotated = ssht_rotate_flm(flm, alpha, beta, gamma);

% Compute inverse harmonic transform.
f_rotated = ssht_inverse(flm_rotated, L, ...
  'Method', args.Method, 'Spin', args.Spin, 'Reality', args.Reality);
