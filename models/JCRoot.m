%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% 	FUNCTION JCRoot.m (MatLab Version 5.1)
% 	THIS FUNCTION CALCULATES ROOTS OF N-TH ORDER JACOBI POLYNOMIAL.  THE EVALUATION
% 	OF THE DERIVATIVES OF THE POLYNOMIALS p IS TAKEN FROM J. VILLADSEN & M. MICHELSEN
%	"SOLUTION OF DIFFERENTIAL EQUATION MODELS BY POLYNOMIAL APPROXIMATION", PRENTICE
%	HALL, NEW JERSEY, 1978.
%
%	Written by	D.D. Do
%					Department of Chemical Engineering
%					University of Queensland, St.Lucia, Qld 4072, Australia
%
%					Date written:	3 January 1998
%-----------------------------------------------------------------------------
%@REF: Villadsen & Michelsen, Section 3.4
function [der, u] = JCRoot(n,n0,n1,al,be)
    % CONSTRAINT ON SOME PARAMETERS
    if n>20; 	disp('The parameter "n" must be smaller than 20'); return; end
    if al<=-1; 	disp('The parameter "al" must be greater than -1'); return; end
    if be<=-1; 	disp('The parameter "be" must be greater than -1'); return; end

    p = jacobiPolynomial(n, al, be);
    %evaluate shifted jacobi polinomial roots
    u = roots(p);							
    if n0==1; 	u = [0 ; u]; 	end
    if n1==1; 	u = [u ; 1]; 	end

    % Sort the roots in ascending order
    u = sort(u);
    %build polynomial with roots u and take derivatives
    der = nodeDerivatesOfPolynomial(u);
    
    global DEBUG
    if DEBUG
         % CONSTRAINT ON SOME PARAMETERS
        if n>20; 	disp('The parameter "n" must be smaller than 20'); return; end
        if al<=-1; 	disp('The parameter "al" must be greater than -1'); return; end
        if be<=-1; 	disp('The parameter "be" must be greater than -1'); return; end

        p1 = jacobiPolynomialMy(n, al, be);
        %evaluate shifted jacobi polinomial roots
        u1 = roots(p1);							
        if n0==1; 	u1 = [0 ; u1]; 	end
        if n1==1; 	u1 = [u1 ; 1]; 	end

        % Sort the roots in ascending order
        u1 = sort(u1);
        %build polynomial with roots u and take derivatives
        der1 = nodeDerivatesOfPolynomialMy(u1);
    end
end

%CALCULATION OF SHIFTED JACOBI POLYNOMIAL COEFICIENTS P
function p = jacobiPolynomial(n, al, be)
    g = ones(1,n+1);
    for i=1:n
        %@REF: Villadsen & Michelsen, eqs. (9) of Chapter 3 
        g(i + 1) = (n-i+1)*(n+i+al+be)*g(i)/(i *(i+be));
    end
    for i=0:n
        %@REF: Villadsen & Michelsen, eqs. (1) of Chapter 3 
        g(i+1) = (-1)^(n-i)*g(i+1);
    end
    p = fliplr(g);
end

%CALCULATION OF DERIVATIVES OF POLYNOMIALS
%@REF: Villadsen & Michelsen, eqs. (51)-(53) of Chapter 3
function der = nodeDerivatesOfPolynomial(polynomialRoots)	
    nt = length(polynomialRoots);
    der = zeros(3, nt);
    for i=1:nt
        der(1,i) = 1; 
        der(2,i) = 0; 
        der(3,i) = 0;
        for j=1:nt;
            if j~=i;
                z = polynomialRoots(i)- polynomialRoots(j);
                der(3,i)= z*der(3,i)+3*der(2,i);
                der(2,i)= z*der(2,i)+2*der(1,i);
                der(1,i)= z*der(1,i);
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%the same as jacobiPolynomial but time consuming
function p = jacobiPolynomialMy(n, al, be)
    syms x;
    jacobi=jacobiP(n, al, be, x);
    syms u;
    jacobi = subs(jacobi, x, 2*u-1);
    w = (1-u)^al * u^be;
    normC = double(int(w*jacobi*jacobi,0,1));
    p = sym2poly(jacobi);
    p = p.*1/sqrt(normC);
end

%The same as nodeDerivatesOfPolynomial but a time consuming.
function der = nodeDerivatesOfPolynomialMy(polynomialRoots)
    polynomial = polynomialWithRoots(polynomialRoots);
    der = nodeDerivates(polynomial, polynomialRoots);
end

%make polynomial with ROOTS as a roots ;)
function p = polynomialWithRoots(roots)
    syms x;
    syms res;
    res=1;
    for i=1:length(roots)
        res = res*(x-roots(i));
    end
    p = sym2poly(res);
end
%calculate 1st,2nd and 3rd derivitieves of polynomial with
%polinomial parameters POLYNOMIAL
function der = nodeDerivates(polynomial, polynomialRoots)
    p = polynomial./polynomial(1);
    p1 = polyder(p);
    p2 = polyder(p1);
    p3 = polyder(p2);
    der = zeros(3, length(polynomialRoots));
    der(1, :) = polyval(p1,polynomialRoots);
    der(2, :) = polyval(p2,polynomialRoots);
    der(3, :) = polyval(p3,polynomialRoots);
end