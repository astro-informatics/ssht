function [ylm] = ssht_ylm(L, theta, phi, varargin)
% ssht_ylm - Compute spherical harmonic functions
%
% Computes spherical harmonic functions using various recursions.
%
% Default usage is given by
%
%   ylm = ssht_ylm(L, theta, phi, <options>)
%
% where L is the harmonic band-limit, theta in [0,pi] is colatitude and phi 
% in [0,2*pi) is longitude.

% Options consist of parameter type and value pairs.  Valid options
% include:
%  'Spin'            = { non-negative integers (default=0) }
%  'Recursion'       = { 'Kostelec'         [3-term recursion, e.g. Kostelec (default)],
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
p.addRequired('theta', @(t) isnumeric(t) && t>=0 && t<=pi);
p.addRequired('phi', @(p) isnumeric(p) && p>=0 && p<=2*pi);
p.addParamValue('Spin', 0, @isnumeric);
p.addParamValue('Recursion', 'Kostelec', @ischar);
p.parse(L, theta, phi, varargin{:});
args = p.Results;

% Compute spherical harmonic functions.








  

%%%%% Matlab

[theta, phi] = s2ea_samples_sphere(B);
[theta_grid, phi_grid] = s2ea_samples_sgrid(B);



lmax = size(flm, 1);


f = zeros(size(theta_grid));
for l=0:lmax-1
    
    Plm_allm = legendre(l,cos(theta));
    c1 = (2*l+1)/(4*pi);
    
    for m=-l:l
        C = sqrt( c1 .* factorial(l-abs(m))/factorial(l+abs(m)) );
        Plm = Plm_allm(abs(m)+1,:);
        [Plm_grid, unused] = meshgrid(Plm, phi);  
        Ylm = C .* Plm_grid .* exp(i.*abs(m).*phi_grid);
        if(m<0) 
            Ylm = (-1).^m * conj(Ylm);  
        end
        f = f + flm(l+1,m+lmax-1+1) .* Ylm;

    end
end


%%%%%%% Numerical Recipes
 
    function s2_ylm_plm(el,m,x) result(plm)
            
      implicit none
      
      integer, intent(in) :: el, m
      real(s2_dp), intent(in) :: x
      real(s2_dp) :: plm

      integer :: i, ll
      real(s2_dp) ::  fact, pll, pmm, pmmp1, somx2
      
      if(m<0 .or. m>el .or. abs(x)>1d0) then
         call s2_error(S2_ERROR_YLM_ARG_INVALID, 's2_ylm_plm')
      end if

      ! Compute Pmm.
      pmm=1.                     
      if(m.gt.0) then
         somx2=sqrt((1.-x)*(1.+x))
         fact=1.
         do i=1,m
            pmm=-pmm*fact*somx2
            fact=fact+2.
         enddo
      endif

      if(el.eq.m) then
         plm=pmm
      else
         ! Compute Pm,m+1.
         pmmp1=x*(2*m+1)*pmm      
         if(el.eq.m+1) then
            plm=pmmp1
         else                    
            ! Compute Pm,l for l > m+1.
            do ll=m+2,el
               pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m)
               pmm=pmmp1
               pmmp1=pll
            enddo
            plm=pll
         endif
      endif
      return

    end function s2_ylm_plm