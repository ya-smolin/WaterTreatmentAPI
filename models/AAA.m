clear
%% data
global al be A B w
global n s Ds Ed Dg
global H M N tau tGrid rGrid
global q Ss x1
%% DONE!: check this result with MATLAB build solver. It works!
al = 0;
be = 2;
x1 = 1;
n = 2.5; %Freundlih isoterm constant
s = 2; %for round
%Problem size
H = 10;  %x
M = 15; %r
N = 20;  %t

tau = 1/(N-1);
tGrid=0:tau:1;

Ds=9.12e-9;
R=0.039;
S0=110;

po=700 * 1000; %my guess
epss=0.7; %my guess
q0=S0^n; %my guess

Dg=po*q0*(1-epss)/(epss*S0);
Ed = Ds*Dg*tau/R^2; %check

q = zeros(H, M, N); %solution
q(:, 1, :) = 0; %(26)
%Ss = @(t)t.^2;
Ss=zeros(x1, N);

% M-1 roots of Jacobi Polinomial
rootsNum = M-2; %-1 for each boundary condition
%% create grid from roots of orthogonal polynomial
yP = YacobiPolynomial(rootsNum, al, be);
rGrid = yP.rootsWithBoundary(1, 1);
lPs = LagrangePolynomials(rGrid);
%% Unit tests section
yP.unitTestOrthogonality();
yP.unitTestYita();
lPs.unitTestPolDer();
lPs.unitTestAB();
Radau.unitTestRDW();
%% Ajm=Lm'(xj), Bjm=Lm''(xj)
A = lPs.derValInNodes1(1);
B = lPs.derValInNodes1(2);
w = Radau.rdwExactInt(rGrid);

tol=1e-10;
iter=1e4;
%'MaxFunctionEvaluations', iter, 'MaxIterations', iter, 'FunctionTolerance', tol, 'StepTolerance', tol,  'OptimalityTolerance', tol
options=optimoptions('fsolve','Display', 'final','Algorithm','levenberg-marquardt');

global rDia tDia rSize tSize eq1Size
rDia=1:M;
tDia = 2:N;
rSize = length(rDia);
tSize = length(tDia);
eq1Size=rSize*tSize;

X = fsolve(@model,zeros(eq1Size, 1), options);
Q = reshape(q(x1, :, :), M, N);
Q(rDia, tDia) = reshape(X, rSize, tSize);
q(x1, rDia, tDia) = Q(rDia, tDia);
sol1=reshape(q(x1, :, :), M, N);
 
[sol2, rr]=getExactSol(N);
plotEvol(sol1, rGrid, sol2, rr);



function Z=model(X)
global al be A B w
global n s Ds Ed Dg
global H M N tau tGrid rGrid rDia tDia rSize tSize eq1Size
global q Ss x1

[q, Ss] = unpack(X);
Z=zeros(rSize, tSize);
Zss=zeros(tSize, 1);
for k = tDia
    for j = rDia
        if j == 1
            eq1 = q(x1, 2, k) - q(x1, 1, k);
        else
            sumR = 0;
            %q(x1, M, k) = Ss(k).^n;
            Zss(k-1) =  q(x1, M, k) - Ss(k).^n;
            for m = 1:M
                sumR = sumR + (s/rGrid(j) * A(j, m)  + B(j, m))*q(x1, m, k);
            end
            eq1 = q(x1, j, k-1) - q(x1, j, k) + tau * Ed * sumR;
        end
        Z(j, k-1) = eq1;
    end
end
Z = reshape(Z, rSize*tSize,1);
Z = [Z;Zss];
end

function [q, Ss] = unpack(X)
global x1 rDia tDia eq1Size rSize tSize
    q(x1, rDia, tDia)=reshape(X(1:eq1Size), rSize, tSize);
    Ss(x1, tDia)=X(eq1Size:eq1Size+tSize);
end

function pack()
end
%grid t have to be the same, for simplicity
function plotEvol(sol1, x1, sol2, x2)
    
    [~, N] = size(sol1);
    tau = 1/(N-1);
    speed = tau; %4sec play
    tInd = @(k)(k-1)*tau;
    ax = axes;
    for k = 1:N
        tic
        title(num2str(tInd(k)));
        if nargin == 4
             plot(ax, x2, sol2(:, k), 'b-');
             hold on;
        end
        plot(ax, x1, sol1(:, k), 'r*');
        axis([0, 1, 0, 1]);
        pause(speed);
        if k ~= N cla; end;
    end
end

function [sol, x]=getExactSol(N)
    global s
    if nargin == 0
        N = 2000;
    end
    x = linspace(0,1,200);
    t = linspace(0,1,N);

    sol = pdepe(s,@(x,t,u,dudx)pdefun(x,t,u,dudx),@icfun,@bcfun,x,t);
    sol = sol';
    
end

function [c,f,s] = pdefun(x,t,u,dudx)
    global Ed
    c=1;
    s=0;
    f = Ed*dudx;
end

function u = icfun(x)
    u = 0;
end

function [pl,ql,pr,qr] = bcfun(xl,ul,xr,ur,t)
global Ed Ss n
ql = 1/Ed;
pl=0;

qr = 0;
pr = ur-Ss(t).^n;
end