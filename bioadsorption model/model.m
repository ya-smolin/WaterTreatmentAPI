function [Z J] = model(X)
global i j;
global h teta;
global C A Psi;
global k n;
global Y kc;
global beta a por ro;

C(i,j+1)=X(1);
A(i,j+1)=X(2);
Psi(i,j+1)=X(3);
Z=zeros(1,3);

[int,jInt]=integrate(i);%границі від 0 до h*(i-1), де i=1..M
Z(1)=C(i,j+1)-exp(-a(teta*(j-1))*h*(i-1))*(Cin(teta*(j-1))+a(teta*(j-1))*int);

r1=(1-por)*beta*teta/(2*ro); r2=teta/(4*Y);
Z(2)=A(i,j+1)-r1*C(i,j+1)+r1*S(i,j+1)+r2*Psi(i,j)*U(i,j+1)...
    +r2*U(i,j)*Psi(i,j+1)+r2*U(i,j+1)*Psi(i,j+1)...
    -r1*(C(i,j)-S(i,j))+r2*U(i,j)*Psi(i,j)-A(i,j);

r3=teta/4; r4=kc*teta/2;
Z(3)= (1-r3*U(i,j)+r4)*Psi(i,j+1)-r3*Psi(i,j+1)*U(i,j+1)...
    -r3*Psi(i,j)*U(i,j+1)+(-1-r3*U(i,j)+r4)*Psi(i,j);

if nargout > 1
    J(1,:)=[1 -a(teta*(j-1))*jInt*k^n*n*A(i,j)^(n-1)   0];
    J(2,1)=-r1;
    J(2,2)=1+r1*dS(i,j+1)+r2*dU(i,j+1)*(Psi(i,j)+Psi(i,j+1));
    J(2,3)=r2*(U(i,j)+U(i,j+1));
    J(3,1)=0;
    J(3,2)=-r3*dU(i,j+1)*(Psi(i,j+1)+Psi(i,j));
    J(3,3)=1-r3*U(i,j)+r4-r3*U(i,j+1);
end
end

function S = S(i,j)
global k A n;
S=(k*A(i,j)).^n;
end

function U = U(i,j)
global u ks;
U=u*S(i,j)/(ks+S(i,j));
end

function dS = dS(i,j)
global k A n;
dS=(k^n)*n*(A(i,j)^(n-1));
end

function dU = dU(i,j)
global u ks;
dU=u*ks*dS(i,j)/(ks+S(i,j))^2;
end

%Бере інтеграл від fint від 0 до (i-1)*h. Якщо i=10, значить інтеграл до
%цього шару по функції fint.
%fint - підінтегральна функція
%upBound = i, по суті верхня межа (i-1)*h, нижня = 0.
function [int, jInt] = integrate(upBound)
global M h;
int=0;
if(ismember(upBound,5:4:M)==1)
    for x=5:4:upBound %Bool integrate O(h^4)
        int=int+(2*h/45)*(7*fint(x-4)+32*fint(x-3)+12*fint(x-2)+32*fint(x-1)+7*fint(x));
    end
    jInt=2*h*7/45;
elseif(ismember(upBound,4:3:M)==1) %Simpson 3/8 integrate O(h^3)
    for x=4:3:upBound
        int=int+(3*h/8)*(fint(x-3)+3*fint(x-2)+3*fint(x-1)+fint(x));
    end
    jInt=3*h/8;
elseif(ismember(upBound,3:2:M)==1) %Simpson integrate O(h^3)
    for x=3:2:upBound
        int=int+(h/3)*(fint(x-2)+4*fint(x-1)+fint(x));
    end
    jInt=h/3;
elseif(upBound==2)%for first - trapetzium integrate O(h^2)
    int=(h/2)*(fint(1)+fint(2));
    jInt=h/2;
elseif(ismember(upBound,8:4:M)==1)%Bool + Simpson 3/8 O(h^3)
    [~,index]=ismember(upBound,8:4:M);
    iter = index;
    x=5;
    while iter>0
        int=int+(2*h/45)*(7*fint(x-4)+32*fint(x-3)+12*fint(x-2)+32*fint(x-1)+7*fint(x));
        x=x+4;
        iter=iter-1;
    end
    int=int+(3*h/8)*(fint(upBound-3)+3*fint(upBound-2)+3*fint(upBound-1)+fint(upBound));
    jInt=3*h/8;
else%6:4:M
    nearest=setdiff(union(4:3:M,3:2:M),5:4:M);
    if(ismember(upBound-2,nearest)==1)
        for x=4:3:upBound-2 %Simson 3/8 + Simpson O(h^3)
            int=int+(3*h/8)*(fint(x-3)+3*fint(x-2)+3*fint(x-1)+fint(x));
        end
        int=int+(h/3)*(fint(upBound-2)+4*fint(upBound-1)+fint(upBound));
        jInt=h/3;
    else %Simson + Simpson 3/8  O(h^3)
        for x=3:2:upBound-3
            int=int+(h/3)*(fint(x-2)+4*fint(x-1)+fint(x));
        end
        int=int+(3*h/8)*(fint(upBound-3)+3*fint(upBound-2)+3*fint(upBound-1)+fint(upBound));
        jInt=3*h/8;
    end
end
end

function y = fint(i)
global A;
global teta j a h n k;
x=h*(i-1);
y=(k*A(i,j+1))^n*exp(a(teta*(j-1))*x);
end
