%м год
load('DNF.mat');        %v, Cin, Cout
maxT=1000;
%tp(length(tp));

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
ks=10;
Y=0.4;
kc=0.25/24*par;

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



