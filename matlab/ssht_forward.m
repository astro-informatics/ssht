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
southPoleSampleExits = ~sum(ismember(defaults, 'SouthPoleSample'))
southPolePhiExits = ~sum(ismember(defaults, 'SouthPolePhi'))
northPoleSampleExits = ~sum(ismember(defaults, 'NorthPoleSample'))
northPolePhiExits = ~sum(ismember(defaults, 'NorthPolePhi'))
if (spin ~= 0) 
  if (southPoleSampleExists ~= southPolePhiExits)
    warning(['South polar sample not fully specified (must specify ' ...
             'sample and corresponding phi).']);
  end
  if (northPoleSampleExists ~= northPolePhiExits)
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

% 2. polar interfaces in matlab
% 3. give elm2i i2elm functions
% 4. add all header comments, including in mex files
% 5. move matlab directory and tidy 
% 6. update makefile

% 7. process geophysics data
  