classdef Radau < handle
    %RADAU Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
      
    end
    methods
         function obj = Radau()
            Radau.unitTestRDW();
        end
    end
    methods (Static)
       
        %
        %TODO: this method works better than original RDW, but have issues
        %      1. for al,be ~= 0 likely works incorrectly
        %      2. it's hard to use because new nodes are generated
        function [w, nodes] = rdwNewNodes(n, n0, n1, al, be)
            if nargin < 3
                error('not enought arguments');
            elseif nargin == 3
                al = 0; 
                be = 0;
            end
            
            N = n+n0+n1;
            %form grid nodes, that is yacobi roots + boundary nodes
            yP = YacobiPolynomial(n, al+n1, be+n0);
            nodes = yP.rootsWithBoundary(n0, n1);

            d = LagrangePolynomials(nodes).pNder();
            w = (2*n + al + be + 1 + n0 + n1)*yP.cn()./(nodes.^n1.*(1-nodes).^n0.*d(1,:)'.^2);
            if n1==1
                w(N)= w(N)/(1+al);
            end
            if n0==1
                w(1) = w(1)/(1+be);
            end
            w = w./sum(w);
        end

        function w = rdwExactInt(nodes, al, be)
            if nargin == 1
                al = 0; 
                be = 0;
            end
            
            lP = LagrangePolynomials(nodes);
            w = zeros(lP.n, 1);
            for i = 1:lP.n
                w(i) = integral(@(x)Radau.weightF(x, al, be).*lP.val(i, x), 0, 1);
            end
        end
        
        function w = weightF(x, al, be)
            w = (1-x).^al .* x.^be;
        end
        
         function pass = unitTestRDW()
            
            n=10;
            n1=1;
            n0=0;
            al=0;
            be=0;
            
            testf = @(x)1 + x.^2 + cos(x).*exp(-x);
            %exact integral using another quadrature
            wF = @(x)Radau.weightF(x, al, be);
            int0 = integral(@(x)wF(x).*testf(x), 0, 1);
            
            [w1, nodes1]=Radau.rdwNewNodes(n, n0, n1, al, be);
            int1 = dot(w1,wF(nodes1).*testf(nodes1));
            
            yP = YacobiPolynomial(n, al+n1, be+n0);
            nodes = yP.rootsWithBoundary(n0, n1);
            w2 = Radau.rdwExactInt(nodes, al, be);
            int2 = dot(w2, testf(nodes));
            
            [answerStr, pass] = equalEpsN(int1, int2, int0);
            disp(['is equal integrals? ', answerStr, ' ' num2str([int1, int2, int0])]);
        end
    end
end

