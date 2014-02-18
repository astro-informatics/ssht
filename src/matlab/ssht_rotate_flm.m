function [flm_rotated] = ssht_rotate_flm(flm, d, alpha, gamma, varargin)
% ssht_rotate_flm - Rotate a function in harmonic space
%
% Rotates a function on the sphere in harmonic space, where
% the rotation is specified by a point (alpha, beta, gamma)
% in SO(3).
%
% Default usage is given by
%
%   [flm_rotated] = ssht_rotate_flm(flm, d, alpha, gamma, <option>)
%
% where flm are the harmonic coefficients of the function to be
% rotated, d are the Wigner small-d functions d_lmn for all el,
% m, n evaluated at beta. They are indexed d(el,m,n), such that
% d has dimensions (L,2*L-1,2*L-1).
% alpha and gamma are the other two rotation angles.
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


% Parse arguments.
p = inputParser;
p.addRequired('flm', @isnumeric);
p.addRequired('d', @isnumeric);
p.addRequired('alpha', @isnumeric);
p.addRequired('gamma', @isnumeric);
p.addParamValue('M', -1, @isnumeric);

p.parse(flm, d, alpha, gamma, varargin{:});

args = p.Results;

L = sqrt(length(flm));

if args.M == -1
    args.M = L;
end

M = args.M;

index = 1;

flm_rotated = zeros(L^2,1);

for el = 0:L-1
    for m = -el:el
        n_max = min(el, M-1);
        for n = -n_max:n_max
            Dlmn = exp(-i*m*alpha) * d(el+1,m+L,n+L) * exp(-i*n*gamma);
            flm_rotated(index) = flm_rotated(index) + ...
                Dlmn * flm(ssht_elm2ind(el,n));
        end
        index = index + 1;
    end
end