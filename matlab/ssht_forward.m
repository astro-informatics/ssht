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
%  'SouthPoleSample' = { double (default=0) }
%  'SouthPolePhi'    = { double (default=0) }
%  'NorthPoleSample' = { double (default=0) }
%  'NorthPolePhi'    = { double (default=0) }
%
% Note that an alternative interface is provided for the MW and MWSS
% methods since these contain samples on the poles (MW sampling
% contains a sample on the South pole, but not North pole; MWSS
% sampling contains a sample on both the North and South poles).
% Instead of passing function values on the poles for all values of
% phi, the poles may be specified by single samples with their
% corresponding phi angle (values for other phi are then related to
% this sample by the rotation of a spin function in its tangent
% plane).
% 
% Author: Jason McEwen (jason.mcewen@epfl.ch)

% Parse arguments.
p = inputParser;
p.addRequired('f', @isnumeric);          
p.addRequired('L', @isnumeric);          
p.addParamValue('Method', 'MW', @ischar);
p.addParamValue('Spin', 0, @isnumeric);
p.addParamValue('Reality', false, @islogical);
p.addParamValue('SouthPoleSample', 0, @isnumeric);
p.addParamValue('SouthPolePhi', 0, @isnumeric);
p.addParamValue('NorthPoleSample', 0, @isnumeric);
p.addParamValue('NorthPolePhi', 0, @isnumeric);
p.parse(f, L, varargin{:});
args = p.Results;

% Check whether individual polar samples specified.
defaults = p.UsingDefaults
southPoleSampleExists = ~sum(ismember(defaults, 'SouthPoleSample'))
southPolePhiExists = ~sum(ismember(defaults, 'SouthPolePhi'))
northPoleSampleExists = ~sum(ismember(defaults, 'NorthPoleSample'))
northPolePhiExists = ~sum(ismember(defaults, 'NorthPolePhi'))
if (args.Spin ~= 0) 
  if (southPoleSampleExists ~= southPolePhiExists)
    warning(['South polar sample not fully specified (must specify ' ...
             'sample and corresponding phi).']);
  end
  if (northPoleSampleExists ~= northPolePhiExists)
    warning(['North polar sample not fully specified (must specify ' ...
             'sample and corresponding phi).']);
  end
end

% Computing forward transform.
flm = ssht_forward_mex(f, L, args.Method, args.Spin, args.Reality, ...
                       (southPoleSampleExists == true), ...
                       args.SouthPoleSample, args.SouthPolePhi, ...
                       (northPoleSampleExists == true), ...
                       args.NorthPoleSample, args.NorthPolePhi);



% 1. polar interfaces in c -- DONE
% 2. polar interfaces in matlab -- DONE
% 3. give elm2i i2elm functions -- DONE
% 4. add all header comments, including in mex files -- DONE

% 4a. create two demos (simple and complex)
% 5. move matlab directory and tidy (create mex directory)
% 6. update makefile

% 7. process geophysics data
  