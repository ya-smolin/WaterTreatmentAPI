%mesuares: m h mg
clear;
global l h M; 
global t teta N;
global C A Psi;
global k n;
global u ks Y kc;
global beta v a por ro;
global C0 psi_0 Cin;

par=1;
parm=1;

l=0.3/parm; h=0.01/parm;
M=round(l/h)+1;
t=100/par; teta=1/par;
N=round(t/teta)+1;

C=zeros(M,N);
A=zeros(M,N);
Psi=zeros(M,N);

k=1/0.0668;
n=1/0.959;

u=0.4*par;
ks=200;
Y=1.35;
kc=0.003*par;

beta=158.4*par;
v=7.2*par/parm;
a=@(t)beta/v;
por=0.6;
ro=352;

psi_0 = 0.0375;
C0=80;
Cin=@(t)C0;

C(:,1) = 0;
A(:,1) = 0;
Psi(:,1) = psi_0;

C(1,:) = C0;
A(1,:) = 1/k*C0^(1/n);
Psi(1,:) = psi_0;




