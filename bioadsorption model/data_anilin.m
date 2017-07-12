%mesuares: m h mg
clear;
global l h M; 
global t teta N;
global C A Psi;
global k n;
global u ks Y kc;
global beta v a por ro;
global C0 psi_0;
global Cin
par=1; %parametrization is turned off
parm=1;

l=0.1/parm; h=0.01/parm;
M=round(l/h)+1;
t=30/par; teta=1/par;
N=round(t/teta)+1;

C=zeros(M,N);
A=zeros(M,N);
Psi=zeros(M,N);

k=1/3.57;%a=(1/k)*S^(1/n)
n=1/0.39;

u=0.25*par;
ks=58.8;
Y=0.6;
kc=0.00003*par;

beta=50*par;
v=5*par/parm;
a=@(t)beta/v;
por=0.6;
ro=300;

psi_0 = 1;
C0=500;
Cin=@(t)C0;

C(:,1) = 0;
A(:,1) = 0;
Psi(:,1) = psi_0;

C(1,:) = C0;
A(1,:) = 1/k*C0^(1/n);
Psi(1,:) = psi_0;





