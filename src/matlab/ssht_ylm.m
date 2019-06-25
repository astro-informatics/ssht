function [ylm] = ssht_ylm(L, theta, phi, varargin)
% ssht_ylm - Compute spherical harmonic functions
%
% Computes spherical harmonic functions for all el and all 0<=|m|<= el  
% using various recursions.
%
% Default usage is given by
%
%   ylm = ssht_ylm(L, theta, phi, <options>)
%
% where L is the harmonic band-limit, theta in [0,pi] is colatitude and phi 
% in [0,2*pi) is longitude.  The results spherical harmonics are indexed
% ylm{ind}(theta,phi), where ind = ssht_elm2ind(el, m).
%
% Options consist of parameter type and value pairs.  Valid options
% include:
%  'Spin'            = { non-negative integers (default=0) }
%  'Recursion'       = { 'Kostelec'         [3-term recursion, e.g. Kostelec (default)],
%                        'Risbo'            [Risbo recursion],
%                        'NumericalRecipes' [Numerical Recipes],
%                        'Matlab'           [Built-in Matlab recursion] }
%
% Author: Jason McEwen (www.jasonmcewen.org)

% SSHT package to perform spin spherical harmonic transforms
% Copyright (C) 2013  Jason McEwen
% See LICENSE.txt for license details

% Parse arguments.
p = inputParser;
p.addRequired('L', @isnumeric);          
p.addRequired('theta', @isnumeric);
p.addRequired('phi', @isnumeric);
p.addParamValue('Spin', 0, @isnumeric);
p.addParamValue('Recursion', 'Kostelec', @ischar);
p.parse(L, theta, phi, varargin{:});
args = p.Results;

% Check size theta and phi identical.
[theta_m, theta_n] = size(theta);
[phi_m, phi_n] = size(phi);
if (theta_m ~= phi_m || theta_n ~= phi_n)
  error('Inconsistent theta and phi data.');
end

% Reshape angles.
theta = theta(:);
phi = phi(:);

% Compute spherical harmonics.
ylm = cell(L^2, 1);
switch args.Recursion
  
  case 'Kostelec'
        
    % Precompute data.
    sqrt_tbl = sqrt([0:2*(L-1)+1])';
    signs = ones(L+1,1);
    signs(2:2:end) = -1;
    
    % Compute spherical harmonics from Wigner functions.
    for itheta = 1:theta_m*theta_n
      dln = zeros(L,1);      
      dlnm1 = zeros(L,1);      
      for el = abs(args.Spin):L-1    
        tmp = dln;
        dln = ssht_dln(dln, dlnm1, ...
          L, el, -args.Spin, theta(itheta), ...
          'SqrtTable', sqrt_tbl, 'SignTable', signs);
        dlnm1 = tmp;
        
        for m = 0:el         
          ind = ssht_elm2ind(el, m);
          ylm{ind}(itheta) = signs(abs(args.Spin)+1) ...
            * sqrt((2*el+1)/(4*pi)) ...
            * dln(m+1) ...
            * exp(1i * m * phi(itheta));
        end
      end
    end       
    
  case 'Risbo'
    
    % Precompute data.
    sqrt_tbl = sqrt([0:2*(L-1)+1])';
    signs = ones(L+1,1);
    signs(2:2:end) = -1;
    
    % Compute spherical harmonics from Wigner functions.
    for itheta = 1:theta_m*theta_n
      dl = zeros(2*L-1,2*L-1);
      for el = abs(args.Spin):L-1
        dl = ssht_dl(dl, L, el, theta(itheta), ...
          'SqrtTable', sqrt_tbl, 'SignTable', signs);  
        for m = 0:el
          ind = ssht_elm2ind(el, m);
          ylm{ind}(itheta) = signs(abs(args.Spin)+1) ...
            * sqrt((2*el+1)/(4*pi)) ...
            * dl(m+L,-args.Spin+L) ...
            * exp(1i * m * phi(itheta));
        end
      end
    end        
    
  case 'NumericalRecipes'
        
    if (args.Spin ~= 0)
      error('Non-zero spin not supported for NumericalRecipes ylm recursion.');
    end
    
    for itheta = 1:theta_m*theta_n
      x = cos(theta(itheta));
    
      for m = 0:L-1
        
        % Compute Pmm.
        pmm = 1.0;        
        somx2 = sqrt((1.-x)*(1.+x));
        fact = 1.0;
        for i = 1:m
          pmm = -pmm * fact * somx2;
          fact = fact + 2.0;
        end        
        ind = ssht_elm2ind(m, m);        
        c1 = (2*m+1)/(4*pi);      
        C = sqrt( c1 .* factorial(m-abs(m))/factorial(m+abs(m)) );
        ylm{ind}(itheta) = C .* pmm .* exp(1i * m * phi(itheta));
        
        % Compute Pm,m+1.
        if (m ~= L-1)
          pmmp1 = x * (2*m+1) * pmm;
          ind = ssht_elm2ind(m+1, m);
          c1 = (2*(m+1)+1)/(4*pi);
          C = sqrt( c1 .* factorial(m+1-abs(m))/factorial(m+1+abs(m)) );
          ylm{ind}(itheta) = C .* pmmp1 .* exp(1i * m * phi(itheta));
        end
        
        % Compute Pm,l for l > m+1.
        for el = m+2:L-1
          pll = (x * (2*el-1) * pmmp1 - (el+m-1) * pmm)/(el-m);
          pmm = pmmp1;
          pmmp1 = pll;
          ind = ssht_elm2ind(el, m);
          c1 = (2*el+1)/(4*pi);      
          C = sqrt( c1 .* factorial(el-abs(m))/factorial(el+abs(m)) );                    
          ylm{ind}(itheta) = C .* pll .* exp(1i * m * phi(itheta));
        end
                        
      end
            
    end
    
  case 'Matlab'
    
    if (args.Spin ~= 0)
      error('Non-zero spin not supported for Matlab ylm recursion.');
    end
    
    for el = 0:L-1      
      Plm = legendre(el, cos(theta));
      c1 = (2*el+1)/(4*pi);      
      for m = 0:el
        ind = ssht_elm2ind(el, m);
        C = sqrt( c1 .* factorial(el-abs(m))/factorial(el+abs(m)) );
        ylm{ind} = C .* Plm(abs(m)+1,:).' .* exp(1i.*abs(m).*phi);        
      end      
    end

  otherwise
    
    error('Invalid recursion method');
    
end

% Reshape output data.
for el = 0:L-1      
  for m = 0:el
    ind = ssht_elm2ind(el, m);
    ylm{ind} = reshape(ylm{ind}, theta_m, theta_n);
  end
end
      