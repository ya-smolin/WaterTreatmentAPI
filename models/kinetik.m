clear
global DEBUG
DEBUG = true;

%%
M = 10;
N = 10;

y = zeros(M, N); %solution
y(1, :) = 1;

% M-1 roots of Jacobi Polinomial
n = M-1;
al = 1;
be = 0;

p = jacobiPolynomial(n, al, be);
u = flip(roots(p));
u = [u; 1];




