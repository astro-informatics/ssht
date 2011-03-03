function [ind] = ssht_elm2ind(el, m)
% ssht_elm2ind - Convert harmonic indices to vector index
%
% Convert (el,m) spherical harmonic indices to ind index to access vector
% of harmonic coefficients (folllowing the Matlab convention, ind is index
% from 1). 
%
% Default usage is given by
%
%   [ind] = ssht_elm2ind(el, m)
%
% Author: Jason McEwen (www.jasonmcewen.org)

% SSHT package to perform spin spherical harmonic transforms
% Copyright (C) 2011  Jason McEwen
% See LICENSE.txt for license details

ind = el .* el + el + m + 1;
