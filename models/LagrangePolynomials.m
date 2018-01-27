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
classdef LagrangePolynomials
    properties
        nodes
        n
    end
    
    methods
        function obj = LagrangePolynomials(nodes)
            obj.nodes = nodes;
            obj.n = length(nodes);
        end
        
        function prod = val(obj, i, x)
            prod = 1;
            for k=1:obj.n
                if(k ~= i)
                    prod = prod .* (x - obj.nodes(k)) ./ (obj.nodes(i) - obj.nodes(k));
                end
            end
        end
        
        %CALCULATION OF DERIVATIVES OF POLYNOMIALS
        %@REF: Villadsen & Michelsen, eqs. (51)-(53) of Chapter 3
        function der = pNder(obj)
            %TODO: lazy init
            der = zeros(3, obj.n);
            for i=1:obj.n
                der(1,i) = 1;
                der(2,i) = 0;
                der(3,i) = 0;
                for j=1:obj.n
                    if j~=i
                        z = obj.nodes(i)- obj.nodes(j);
                        der(3,i)= z*der(3,i)+3*der(2,i);
                        der(2,i)= z*der(2,i)+2*der(1,i);
                        der(1,i)= z*der(1,i);
                    end
                end
            end
        end
        
        %calculate 1st,2nd and 3rd derivitieves of polynomial
        function der = pNderSlow(obj)
            p = poly(obj.nodes);
            p1 = polyder(p);
            p2 = polyder(p1);
            p3 = polyder(p2);
            der = zeros(3, obj.n);
            der(1, :) = polyval(p1,obj.nodes);
            der(2, :) = polyval(p2,obj.nodes);
            der(3, :) = polyval(p3,obj.nodes);
        end
        
        function pass = unitTestPolDer(obj)
            pass = true;
            der1 = obj.pNder();
            der2 = obj.pNderSlow();
            if sum(sum(abs(der1-der2))) > 50
                pass = false;
                warning('jacobi polinomial node derivatives can be incorrect');
            end
        end
        
        % Aij=Lj'(xi), Bij=Lj''(xi)
        function v = derValInNodes1(obj, id)
            d = pNder(obj);
            % DIMENSIONALIZE THE OUTPUT VECTOR
            v = zeros(obj.n);
            % START THE CALCULATION
            for i = 1:obj.n
                for j = 1:obj.n
                    if j == i
                        if id==1
                            v(i,i) = d(2,i)/d(1,i)/2;								% Eq.(40),Chapter 3, V&M
                        else
                            v(i,i) = d(3,i)/d(1,i)/3;								% Eq.(41),Chapter 3, V&M
                        end
                    else
                        v(i,j)= d(1,i)/d(1,j)/(obj.nodes(i)-obj.nodes(j));				% Eq.(42),Chapter 3, V&M
                        if id == 2
                            v(i,j) = v(i,j)*(d(2,i)/d(1,i) - 2/(obj.nodes(i)-obj.nodes(j)));	% Eq.(43),Chapter 3, V&M
                        end
                    end
                end
            end
            
        end
        
        %A = L(:,:,1); B = L(:,:,2); Aij=Lm'(xi), Aij=Li'(xj), Bij=Lm''(xi) Bij=Li''(xj)
        function L = derValInNodes2(obj)
            der = pNder(obj);
            L = zeros(obj.n, obj.n, 2);
            for m=1:obj.n
                for i=1:obj.n
                    if i == m
                        L(m, i, 1) = 1/2 * der(2, i)/ der(1,i);
                        L(m, i, 2) = 1/3 * der(3,i)/ der(1,i);
                    else
                        L(m, i, 1) = (der(1,i)/der(1,m) -  lmInd(m, i))./(obj.nodes(i)-obj.nodes(m));
                        L(m, i, 2) = (der(2,i)/der(1,m) -  2*L(m, i, 1))./(obj.nodes(i)-obj.nodes(m));
                    end
                end
            end
            
            function prod = lmInd(m, i)
                if i == m
                    prod = 1;
                else
                    prod = 0;
                end
            end
        end
        
        function pass = unitTestAB(obj)
            pass = true;
            A1 = derValInNodes1(obj, 1);
            B1 = derValInNodes1(obj, 2);
            
            L = derValInNodes2(obj);
            A2 = L(:,:,1)';
            B2 = L(:,:,2)';
            
            if sum(sum(abs(A1-A2+B1-B2))) > 0.1
                pass = false;
                warning('something wrong :(');
            end
        end
        
        function w = rdw(obj,al,be)
            n0 = (obj.nodes(1) == 0);
            n1 = (obj.nodes(end) == 1);
            N = obj.n - n0 - n1; %interior points
            
            yP = YacobiPolynomial(N, al+n1, be+n0);
            nnodes = [yP.u; 1];
            nLP = LagrangePolynomials(nnodes);
            d = nLP.pNder();
            x = nLP.nodes;
            w = (2*N + al + be + 1 + n0 + n1)*yP.cn()./(x.^n1.*(1-x).^n0.*d(1,:)'.^2);
            if n1==1
                w(obj.n)= w(obj.n)/(1+al);
            end
            if n0==1
                w(1) = w(1)/(1+be);
            end
            w = w./sum(w);
        end
        
        function w = rdwMy(obj, al, be)
            w = zeros(obj.n,1);
            wf = @(x)(1-x).^al .* x.^be;
            for i = 1:obj.n
                w(i) = integral(@(x)wf(x).*val(obj, i, x), 0, 1);
            end
        end
        
        function pass = unitTestRDW(obj)
            pass = true;
            testf = @(x)x.^2;
            testF01 = 1/3;
            
            w1=rdw(obj,0,0);
            int1 = w1'*testf(obj.nodes);
            
            w2=rdw00slow(obj);
            int2 = w2'*testf(obj.nodes);
            
            if abs(int2 - testF01)  > 0.001 || abs(int1 - testF01)  > 0.001
                pass = false;
                warning(['incorrect quadrature calc ' num2str(int1) ' ' num2str(int2) ' ' num2str(testF01)]);
            end
            if sum(abs(w1-w2)) > 0.1
                pass = false;
                warning('different methods different result in weights calculation');
            end
        end
    end
end
