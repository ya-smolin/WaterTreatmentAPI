clear;
global i j;
global l h M;
global t teta N;
global C A Psi;
global k n;
global u ks Y kc;
global beta a por ro;
global tp vp CinM CoutM;
%м год
%=====data begin===========================================================
load('DNF.mat'); %for global functions: v, Cin, Cout
maxT=1500;%tp(length(tp));

par=10^6;%на 10^6 год
parm=10;
l=0.085/parm; h=0.01/parm;
M=round(l/h)+1;
t=maxT/par; teta=10/par;
N=round(t/teta)+1;

C=zeros(M,N);
A=zeros(M,N);
Psi=zeros(M,N);

k=1/exp(4.89);
n=1/0.2922;

% u=1/24*par;
% ks=10;
% Y=0.4;
% kc=0.25/24*par;

beta=35*par;
a=@(t)beta/(v(t)*par/parm);
por=0.6;
ro=520;

psi_0 = 0;

C(:,1) = 0;
A(:,1) = 0;
Psi(:,1) = psi_0;

for j=1:N
    C(1,j) = Cin(teta*(j-1));
    A(1,j) = 1/k*Cin(teta*(j-1))^(1/n);
end
Psi(1,:) = psi_0;
%===data end===============================================================


fit=zeros(length(u),length(ks),length(Y),length(kc));
fitr2=fit;

CexitM=zeros(length(u),length(ks),length(Y),length(kc),N);

yT=zeros(1,N);
x=0:teta*par:maxT;
for j=1:N
    yT(j)=Cout(x(j));
end

kci=1;
kc=0.01*par:0.5*par:1.01*par;
Y=0.1:0.4:0.9;
for kcp=kc
    
    
    Yi=1;
    for Yp=Y
      AAA{kci,Yi}=[kcp;Yp];
    Yi=Yi+1;
    end
kci=kci+1;
end
