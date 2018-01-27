%%Ravindran V, Badriyha BN, Pirbazari M, KimSH.
%Modeling of bioactive carbon adsorbers: a hybrid weighted residual-finite difference numerical technique.
%Appl Math Comp 1996;76(2–3):99–131.
global DEBUG
DEBUG = 0;
clear
%% Kinetic equation (adsorbent phase material balance) (20)
%dq(x,r,t)/dt = Ed/r^2*d/dr[r^2*dq(x,r,t)/dt]
%q(x,r,0)=0; q(x, 0, t) = 0; q(x, 1, t) = S(x,t)^n (26)

%Let's calculate just for x==0
%S(0,t) = 1; (25)
Sxt = 1; % <-- link to the next equation, temporal hardcode
n = 2.5; %Freundlih isoterm constant

%Problem size
H = 3;  %x
M = 10; %r
N = 5;  %t

q = zeros(H, M, N); %solution
q(:, 1, :) = 0; %(26)
q(:, M, :) = Sxt.^n;

% M-1 roots of Jacobi Polinomial
rootsNum = M-2; %-1 for each boundary condition
al = 1;
be = 0;

yP = YacobiPolynomial(rootsNum, al, be);
disp(['check Orthogonality is ' num2str(yP.unitTestOrthogonality())]);
nodes = [0; yP.u; 1]; %add boundary conditions

lPs = LagrangePolynomials(nodes);
disp(['check node derivatives ' num2str(lPs.unitTestPolDer())]);
disp(['check AB matrixes ' num2str(lPs.unitTestAB())]);
%% Aij=Lj'(xi), Bij=Lj''(xi)
%A = lPs.derValInNodes1(1);
%B = lPs.derValInNodes1(2);

lPs.unitTestRDW()
