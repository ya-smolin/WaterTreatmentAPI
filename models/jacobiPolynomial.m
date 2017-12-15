%%
% Shifted jacobi polynomial coeficients,
% x=2*u-1, u ? (0,1), x ? (-1,1)
function p = jacobiPolynomial(N, al, be)
    global DEBUG
    if DEBUG
      SAMPLE_NUM = 5;
      SAFE_JACOBI_RANGE_MAX = 20;
      SAFE_JACOBI_RANGE_MIN = 3;
      ij = randi([SAFE_JACOBI_RANGE_MIN SAFE_JACOBI_RANGE_MAX], 2, SAMPLE_NUM);
      for sample = 1:SAMPLE_NUM
          i = ij(1, SAMPLE_NUM);
          j = ij(2, SAMPLE_NUM);
          if i == j 
              continue; 
          end
          intVal = checkOrthogonality(i, j, al, be);
          if(intVal > eps('double') || intVal < -eps('double'))
              warning('incorrect jacobi polynomials')
          end
      end
    end
    p = jacobiPolynomialUnit(N, al, be);
    
end
%%
function p = jacobiPolynomialUnit(N, al, be)
    if N > 23
        warning('Can be incorrect after 23th order');
    end
    p = zeros(N+1,1);
    p(1) = 1;
    for i = 1:N
        p(i+1) = -p(i)*(N-i+1)*(N+i+al+be)/(i*(i+be));
    end
    p = flip(p);
end
%%
function intVal = checkOrthogonality(i, j, al, be)
    syms x;
    w = (1-x)^al * x^be;
    f = poly2sym(jacobiPolynomial(i, al, be),x);
    g = poly2sym(jacobiPolynomial(j, al, be),x);
    intVal = double(int(w*f*g,0,1));
end