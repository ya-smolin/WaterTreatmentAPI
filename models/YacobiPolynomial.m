
classdef YacobiPolynomial < handle
    % Orthogonal Polinomial respectively to the weight function W(x)=x^be(1-x)^al
    % in the interval [0,1]. In literature they called SHIFTED due to the
    % transformation x=2*u-1 on the domains x[-1, 1] to domains u[0, 1].
    properties
        al, be %weight function  W(x)=x^be(1-x)^al params, -1 < al, be
        n      %polynomial degree == roots num
        p      %polynomial params p'*[1, x, x^2, ..., x^n]
        u      %polynomial roots in interval [0, 1]
    end
    
    properties (GetAccess = private)
        isPcalculate = false;
    end
    
    methods
        
        function obj = YacobiPolynomial(n, al, be)
            if n > 23
                warning('Can be incorrect after 23th order');
            end
            if al < -1 || be < -1
                warning('incorrect weight function parameters');
            end
            obj.al = al;
            obj.be = be;
            obj.n = n;
            obj.p = polynomialParams(obj);
            obj.u = polynomialRoots(obj);
        end
        
        %CALCULATION OF SHIFTED JACOBI POLYNOMIAL COEFICIENTS P
        function p = polynomialParams(obj)
            g = ones(1,obj.n+1);
            for i=1:obj.n
                %@REF: Villadsen & Michelsen, eqs. (9) of Chapter 3
                g(i + 1) = (obj.n-i+1)*(obj.n+i+obj.al+obj.be)*g(i)/(i *(i+obj.be));
            end
            for i=0:obj.n
                %@REF: Villadsen & Michelsen, eqs. (1) of Chapter 3
                g(i+1) = (-1)^(obj.n-i)*g(i+1);
            end
            p = fliplr(g);
        end
        
        function u = polynomialRoots(obj)
            u = sort(roots(obj.p));
            if sum(abs(u) > 1)~=0 || sum(u < 0) ~= 0
                warning(['roots is out of (0,1): ' num2str(u')]);
            end
        end
        
        function pass = unitTestOrthogonality(obj)
            warning('off', 'MATLAB:integral:MaxIntervalCountReached');
            ansStr = 'Yes';
            pass = true;
            SAMPLE_NUM = 5;
            SAFE_JACOBI_RANGE_MAX = 20;
            SAFE_JACOBI_RANGE_MIN = 3;
            ij = randi([SAFE_JACOBI_RANGE_MIN SAFE_JACOBI_RANGE_MAX], 2, SAMPLE_NUM);
            for sample = 1:SAMPLE_NUM
                i = ij(1, SAMPLE_NUM);
                j = ij(2, SAMPLE_NUM);
                if i ~= j
                    intVal = checkOrthogonality(obj, i, j);
                    [ansStr, pass] = equalEps(intVal, 1e-10);
                    break;
                end
            end
            disp(['your Jacobi polynomials are orthogonal? :', ansStr]);
        end
        
        function intVal = checkOrthogonality(obj, i, j)
            yP1 = YacobiPolynomial(i, obj.al, obj.be);
            yP2 = YacobiPolynomial(j, obj.al, obj.be);
            intVal = integral(@(x)weightF(obj, x).*yP1.val(x).*yP2.val(x),0,1);
            %w = warning('query','last');
            %id = w.identifier;
            %warning('off',id);
        end
        
        function w = weightF(obj, x)
            w = (1-x).^obj.al .* x.^obj.be;
        end
        
        function v = val(obj, x)
            pF = obj.p;
            v = polyval(pF, x);
        end
        
        function yita1 = yNdirect(obj, i)
            N = obj.n;
            yita1 = nchoosek(N,i)*gamma(N+i+obj.al+obj.be+1)*gamma(obj.be+1)...
                /(gamma(N+obj.al+obj.be+1)*gamma(i+obj.be+1));
        end
        
        function yita2 = yNrecursive(obj, k)
            N = obj.n;
            y = zeros(k+1,1);
            y(1) = 1;
            for i = 1:k
                y(i+1) =(N-i+1)*(N+i+obj.al+obj.be)/(i*(i+obj.be))*y(i);
            end
            yita2 = y(k+1);
        end
        
        function yita3 = yNparams(obj, i)
            N = obj.n;
            yita3 = obj.p(N-i+1)*(-1)^(N-i);
        end
        
        function ret = cn(obj)
            yita = yNparams(obj, obj.n);
            ret = Cn(obj)./yita.^2;
        end
        
        function ret = Cn(obj)
            N = obj.n;
            ret = gamma(obj.be+1).^2*factorial(N)*gamma(N + obj.al + 1)...
                /(gamma(N+obj.be+1)*gamma(N+obj.al+obj.be+1)*(2*N+obj.al+obj.be+1));
        end
        
        function unitTestYita(obj)
            i = obj.n;
            yita1 = yNdirect(obj, i);
            yita2 = yNrecursive(obj, i);
            yita3 = yNparams(obj, i);
            
            disp(['are yNdirect, yNrecursive and yNparams the same? ', equalEpsN(yita1, yita2, yita3)...
                ' ' num2str([yita1, yita2, yita3])]);
        end
        
        function nodes = rootsWithBoundary(obj, n0, n1)
            nodes = obj.u;
            if n1==1
                nodes = [nodes; 1];
            end
            if n0==1
                nodes = [0; nodes];
            end
        end
    end
    
    
    
end

