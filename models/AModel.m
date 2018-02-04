%%Ravindran V, Badriyha BN, Pirbazari M, KimSH.
%Modeling of bioactive carbon adsorbers: a hybrid weighted residual-finite difference numerical technique.
%Appl Math Comp 1996;76(2–3):99–131.
global DEBUG
DEBUG = 0;
clear
% xi=1;
% Sxi = @(t)1;
% qxi = kineticEq(xi, Sxi);

testf = @(x)x.^2 ;
al = 0;
be = 2;
%exact integral using another quadrature
wF = @(x)Radau.weightF(x, al, be);
int0 = integral(@(x)wF(x).*testf(x), 0, 1);
            
yP = YacobiPolynomial(10, al, be);
rGrid = yP.rootsWithBoundary(1, 1);
w2 = Radau.rdwExactInt(rGrid, al, be);
int2 = dot(w2, testf(rGrid));     
[answerStr, pass] = equalEpsN(int2, int0);
disp(['is equal integrals? ', answerStr, ' ' num2str([int2, int0])]);
%% Kinetic equation (adsorbent phase material balance) (20)
%dq(x,r,t)/dt = Ed/r^s*d/dr[r^s*dq(x,r,t)/dt]
%dq(x,r,t)/dt = (s/r)*dq(x,r,t)/dr + Ed*ddq(x,r,t)/dr^2
%
%q(x,r,0)=0; dq(x, 0, t)/dr = 0; q(x, 1, t) = S(x,t)^n (26)
function qxi=kineticEq(xi, Sxi)
%Let's calculate just for x==0
%S(0,t) = 1; (25)

n = 2.5; %Freundlih isoterm constant
s = 2; %for round
%Problem size
H = 3;  %x
M = 18; %r
N = 300;  %t

tau = 1/N;
tGrid=tau:tau:1;

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
q(:, M, :) = Sxi(tGrid).^n;

% M-1 roots of Jacobi Polinomial
rootsNum = M-2; %-1 for each boundary condition
al = 1;
be = 0;

yP = YacobiPolynomial(rootsNum, al, be);

rGrid = yP.rootsWithBoundary(1, 1);
lPs = LagrangePolynomials(rGrid);

%% Unit tests section
yP.unitTestOrthogonality()
yP.unitTestYita();
lPs.unitTestPolDer();
lPs.unitTestAB();
Radau.unitTestRDW();
%% Ajm=Lm'(xj), Bjm=Lm''(xj)
A = lPs.derValInNodes1(1);
B = lPs.derValInNodes1(2);
w = Radau.rdwExactInt(rGrid);

xi = 1;% <- hardcode, connection with another equation
Alocal = zeros(M-1, M-1);
b = zeros(M-1, 1);
for k = 2:N
    for j = 2:M-1
        b(1) = 0;
        Alocal(1, 1:2) = [-1 1];
        b(j) = (1/tau)*q(xi, j, k-1) - (s/rGrid(j)*A(j,M)-Ed*B(j, M))*Sxi(j)^n;
        for m = 1:M-1
            Alocal(j, m) = (s/rGrid(j)*A(j,m)-Ed*B(j, m));
            if m == j
                Alocal(j, m) =  Alocal(j, m) + 1/tau;
            end
        end
    end
    debugSol = Alocal \ b;
    q(xi, 1:M-1, k) = debugSol;
end
qxi = q(xi, :, :);
qxi = reshape(qxi, M, N);
plotQ(q, rGrid, xi);
end



function plotQ(q, rGrid, xi)
    [~, ~, N] = size(q);
    tau = 1/N;
    tInd = @(k)(k-1)*tau;
    ax = axes;
    for k = 1:round(N/5)
        title(num2str(tInd(k)));
        plot(ax, rGrid, q(xi, :, k), 'r-');
        axis([0, 1, 0, 1]);
        hold on;
        pause(0.05);
    cla;
    end
    close all;
end