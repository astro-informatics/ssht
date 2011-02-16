function [f, varargout] = ssht_inverse(flm, L, varargin)
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
% Note that an alternative interface is provided for the MW and MWSS
% methods since these contain samples on the poles (MW sampling
% contains a sample on the South pole, but not North pole; MWSS
% sampling contains a sample on both the North and South poles).
% Instead of requesting function values on the poles for all values of
% phi, the poles may be specified by single samples with their
% corresponding phi angle (values for other phi are then related to
% this sample by the rotation of a spin function in its tangent
% plane).  For the MW method this interface is called through the
% usage
%
%   [f, f_sp, phi_sp] = ssht_inverse(flm, L, <options>)
%
% where f does not contain samples on the South pole, f_sp is the
% South pole sample and phi_sp is its corresponding phi value.  For
% the MWSS method this interface is called through the usage
%
%   [f, f_sp, phi_sp, f_np, phi_np] = ssht_inverse(flm, L, <options>)
%
% where f_np is the North pole sample and phi_np is its corresponding
% phi value.
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
southPoleExists = (nargout >= 3);
northPoleExists = (nargout >= 5);
args = p.Results

% Computing inverse transform.
[f, f_sp, phi_sp, f_np, phi_np] = ...
    ssht_inverse_mex(flm, L, args.Method, args.Spin, args.Reality, ...
                     southPoleExists, northPoleExists);
                  
% Set optional outputs.
if ~(nargout == 1 || nargout == 3 || nargout == 5)
  error('Invalid number of output arguments.');
end
if nargout >= 3
  varargout(1) = {f_sp};
  varargout(2) = {phi_sp};
end
if nargout >= 5
  varargout(3) = {f_np};
  varargout(4) = {phi_np};
end



  