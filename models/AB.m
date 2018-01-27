%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% 	FUNCTION AB.M (MatLab Version 5.1)
%	THIS FUNCTION CALCULATES THE FIRST & SECOND DERIVATIVE MATRICES FOR
%	THE ORTHOGONAL COLLOCATION METHOD
%					The algorithm is taken from Villadsen & Michelsen (1978) (V&M)
%					"Solution of Differential Equation Models by Polynomial Approximation"
%					Prentice Hall, New Jersey, 1978, eqs. (40) to (43) pgs 120, 121.
%
%					Date written: 15 October 1997
%
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Nomenclature:
%	d			:	dif matrix from jcob function program
%	id			:	index on the nature of output
%					=1 for calculation of the first derivative matrix
%					=2 for calculation of the second derivative matrix
%	n			:	number of interior collocation points
%	n0			:	=1 if x=0 is included, otherwise =0
%	n1			:	=1 if x=1 is included, otherwise =1
%	u			:	interpolation points of length (n+n0+n1)
%
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% USAGE:
%	v 	= AB(n,n0,n1,id,d,u)
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++	
 function v = AB(n,n0,n1,id,d,u)
	nt		= n+n0+n1;

% DIMENSIONALIZE THE OUTPUT VECTOR
	v		= zeros(nt);

% START THE CALCULATION 
	for i=1:nt;
		for j=1:nt;
			if j==i;
				if id==1;
					v(i,i)	= d(2,i)/d(1,i)/2;								% Eq.(40),Chapter 3, V&M
                else
					v(i,i)	= d(3,i)/d(1,i)/3;								% Eq.(41),Chapter 3, V&M
				end;
			else
				v(i,j)		= d(1,i)/d(1,j)/( u(i)-u(j) );				% Eq.(42),Chapter 3, V&M
				if id==2; 
					v(i,j) = v(i,j)*( d(2,i)/d(1,i)-2/(u(i)-u(j)) );	% Eq.(43),Chapter 3, V&M
				end;
			end
		end
	end;
    
 end   

function v = ABmy(n,n0,n1,id,d,u)
	nt		= n+n0+n1;
    P = zeros(nt, 3);
    P = d';
    xV = u;
    L = lmIndex(xV, P);
     
    if id == 1
       v =  L(:,:,1);
    else
       v =  L(:,:,2);
    end
end 
        
function L = lmIndex(xV, P)
    M = length(xV);
    L = zeros(M, M, 2);
    for m=1:M
        for i=1:M
            if i == m
                L(i, m, 1) = 1/2 * P(i,2)/ P(i,1);
                L(i, m, 2) = 1/3 * P(i,3)/ P(i,1);
            else
                L(i, m, 1) = (P(i,1)/P(m,1) -  lmInd(i, m))./(xV(i)-xV(m));
                L(i, m, 2) = (P(i,2)/P(m,2) -  2*L(i, m, 1))./(xV(i)-xV(m));
            end
            
        end
    end
end

function prod = lmInd(i, m)
    if i == m
        prod = 1;
    else
        prod = 0;
    end
end
