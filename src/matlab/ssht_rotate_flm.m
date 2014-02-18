function [flm_rotated] = ssht_rotate_flm(flm, dlm, alpha, gamma, varargin)
% ssht_rotate_flm - Rotate a function in harmonic space
%
% Rotates a function on the sphere in harmonic space, where
% the rotation is specified by a point (alpha, beta, gamma)
% in SO(3).
%
% Default usage is given by
%
%   [flm_rotated] = ssht_rotate_flm(flm, dlm, alpha, gamma, <option>)
%
% where flm are the harmonic coefficients of the function to be
% rotated, dlm is the Wigner small-d functions d_lm evaluated
% at beta and alpha and gamma are the other two rotation angles.
% flm_rotated are the harmonic coefficients of the resulting
% function after rotation.
% The band-limit L is implied by the size of the flm, which
% should be L*L.
%
% Options consist of parameter type and value pairs.
% Valid options include:
%
%  'M' = Non-negative integer, not greater than L
%        (defaults to L). Azimuthal band-limit.
%
% Authors: Martin Büttner (m.buettner.d@gmail.com)
%          Jason McEwen (www.jasonmcewen.org)

% SSHT package to perform spin spherical harmonic transforms
% Copyright (C) 2011-2014 Jason McEwen
% See LICENSE.txt for license details
