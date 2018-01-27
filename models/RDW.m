%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Function RDW MatLab version 5.1
% THIS FUNCTION CALCULATES THE NORMALISED RADAU WEIGHTS
% Ref: J. Villadsen & M. Michelsen, "Solution of Differential Equation Models
%			by Polynomial Approximation", Prentice Hall, New Jersey, 1978
%			Chapter 3, eqs. (83)-(85)
%
% Written by 	D. D. Do
%					Department of Chemical Engineering
%					University of Queensland
%					St. Lucia, Qld 4072, Australia
%
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% 	Nomenclature:
%	al			:	exponent of the integrand factor (1-x)
%	be			:	exponent of the integrand factor x
%	d			:	dif matrix from JCRoot function routine
%	n			:	number of interior collocation point
%	n0			:	=1 if x=0 is included, otherwise =0
%	n1			:	=1 if x=1 is included, otherwise =0
%	u			:	interpolation points, having length (n+n0+n1)
%
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% USAGE:
%	w = RDW(n,n0,n1,al,be,u,d)
%
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function w = RDW(n,n0,n1,al,be,u,d)
   nt		= n+n0+n1;
	w		= zeros(1,nt);
%
	if n0==0 & n1==1;
      w		= 1./u./d(1,:)'.^2;
      w(nt)	= w(nt)/(1+al);
   end
	if n0==1 & n1==0;
      w		= 1./(1-u)./d(1,:)'.^2;
      w(1)	= w(1)/(1+be);
   end
	if n0==1 & n1==1;
      w		= 1./d(1,:)'.^2;
      w(1)	= w(1)/(1+be);
      w(nt)	= w(nt)/(1+al);
   end
%
	w	= w./sum(w);
