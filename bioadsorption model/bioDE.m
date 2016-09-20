function CC=bioDE(u_,ks_,Y_,kc_)
%=====dataspace begin======================================================
global i j;
global l h M;
global t teta N;
global C A Psi;
global k n;
global u ks Y kc;
global beta a por ro;
global tp vp CinM CoutM;
%м год
%=====dataspace end========================================================
u = u_;
ks = ks_;
Y = Y_;
kc = kc_;
%options=optimset('Jacobian','on','MaxFunEvals',5000,'MaxIter', 5000,'TolX',1e-12,'TolFun',1e-12);
options=optimset('Jacobian','on');
for j = 1:N-1 %для всього часу
    for i = 2:M %для всього шару вугілля при і=1 - це граничні умові задані раніше
        if (j == 1)
            X = fsolve(@model,[C(1,1);A(1,1);Psi(1,1)],options); %для першої години - наближення це граничні умови
        else
            X = fsolve(@model,[C(i,j);A(i,j);Psi(i,j)],options); %для наступних - попередній шар
        end
        %перенесення розв'язку в матрицю вихідних данних
        if(isreal(X(1))&&X(1)>0)
            C(i,j+1)=X(1);
        else
            C(i,j+1)=0;
        end
        
        if(isreal(X(2))&&X(2)>0)
            A(i,j+1)=X(2);
        else
            A(i,j+1)=0;
        end
        
        if(isreal(X(3))&&X(3)>0)
            Psi(i,j+1)=X(3);
        else
            Psi(i,j+1)=0;
        end
        
    end
end
CC=C(M,:);
end