clear
global DEBUG
DEBUG = true;

%Grid size
M = 10;
N = 10;

y = zeros(M, N); %solution
y(:, 1) = 0;
y(M, :) = 1;


% M-1 roots of Jacobi Polinomial
n = M-1;
al = 1;
be = 0;

p = jacobiPolynomial(n, al, be);
xV = flip(roots(p));
xV = [xV; 1];
Tmax = 1;
theta = Tmax / N;
t = 0:theta:Tmax;

P = nodeDerivatesOfPolynomial(xV);
L = lmIndex(xV, P');
A = L(:,:, 1);
B = L(:,:, 2);

w = RDW(n,0,1,al-1,be,xV,P);
toc
funTest = @(x)x.^2;
int1 = w'*funTest(xV);


ww = zeros(M,1);
for i = 1:M
    ww(i) = integral(@(x)lmUnit(x, i, xV), 0, 1);
end
int2 = ww'*funTest(xV);
disp(['int ' num2str([int1 int2 1/3])])

% n=9;
% n0=1;
% n1=1;
% A				= AB(n,n0,n1,1,P,xV)';
% B				= AB(n,n0,n1,2,P,xV)';
    
    
yL = zeros(M-1, M-1);
b = zeros(M-1,1);
%for each layer in time
for j = 2:N
    yL = zeros(M-1, M-1);
    for i = 1:M-1
        for m =1:M-1
            yL(i, m) = 4*xV(i)*B(m,i)+4*A(m,i);
        end
        
        alha = 1 / theta;
        yL(i, i) = yL(i, i) - alha;
        b(i) = - alha*y(i,j-1) - (4*xV(i)*B(M,i)+4*A(M,i))*y(M,j);
    end
    y(1:M-1,j) = yL \ b;
end

%Runge Kutta
tspan = [0, 1];
y0 = zeros(M-1, 1);
yout = ode15s(@(t,y)odefun(t, y, xV, A, B), tspan, y0);
yy = yout.y;

function ans = odefun(t, y, x, A, B)
    y = [y; 1];
    M = length(x);
    ans = zeros(M-1, 1);
    for i = 1:M-1
        ans(i) = sum((4*x(i).*B(:,i)+4*A(:,i)).*y);
    end
end

function plotAll()
    if DEBUG
        rootSize = length(xV);
        baloneyX = zeros(rootSize, 1);
        baloneyY = zeros(N+1, 1);
        plot(baloneyX, xV, 'r*');
        hold on;
        plot(baloneyX+theta, xV, 'ro');
        %plot(t, baloneyY, 'ro');
        ylabel("u");
        xlabel("t");
        set(gca,'XAxisLocation', 'top');
        set(gca,'Ydir','reverse')
        set(gca,'Xdir','normal')
        set(gca,'ytick', xV)
        grid on;

        text(0,xV(end), '  M');
        text(t(end), 0, '  N');

        randomi = floor(M/2 -2);
        randomj = 2;
        text(0,xV(randomi), '  I');
        text(t(randomj), 0, '  J');
        plot(t(randomj), xV(randomi), 'rp')

        plot(t, 1, 'r*')
    end
end
        


function prod = lmInd(m, i)
    if i == m
        prod = 1;
    else
        prod = 0;
    end
end

function prod = lmUnit(x, m, xV)
prod = 1;
M = length(xV);
    for k=1:M
        if(k ~= m)
            prod = prod .*  (x -  xV(k)) ./ (xV(m) - xV(k));
        end
    end
end
