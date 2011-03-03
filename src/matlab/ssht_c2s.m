function [theta, phi, r] = ssht_c2s(x, y, z)
% ssht_c2s - Convert cartesian to spherical coordinates
%
% Convert cartesian to spherical coordinates.  
%
% Default usage is given by
%
%   [theta, phi, r] = ssht_c2s(x, y, z)
%
% Author: Jason McEwen (www.jasonmcewen.org)

% SSHT package to perform spin spherical harmonic transforms
% Copyright (C) 2011  Jason McEwen
% See LICENSE.txt for license details

TOL = 1e-10;

r = sqrt(x.^2 + y.^2 + z.^2);

% Initialise angles.
theta = 0;
phi = 0;

% If radius zero, then angles arbitrary.  Leave set to zero.
% Otherwise ...
if abs(r - TOL) > 0

    theta = atan2(sqrt(x.^2 + y.^2), z);

    % If theta zero, then phi arbitrary.  Leave set to zero.
    % Otherwise ...
    if abs(theta - TOL) > 0
        phi = atan2(y, x);
        phi = mod(phi, 2.*pi);
    end
    
end
