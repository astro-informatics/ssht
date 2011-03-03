function [x, y, z] = ssht_s2c(theta, phi, r)
% ssht_s2c - Convert spherical to cartesian coordinates
%
% Convert cartesian to spherical coordinates.  
%
% Default usage is given by
%
%   [x, y, z] = ssht_s2c(theta, phi, r)
%
% Author: Jason McEwen (www.jasonmcewen.org)

% SSHT package to perform spin spherical harmonic transforms
% Copyright (C) 2011  Jason McEwen
% See LICENSE.txt for license details

x = r .* sin(theta) .* cos(phi);
y = r .* sin(theta) .* sin(phi);
z = r .* cos(theta);