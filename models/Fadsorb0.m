function yprime=f(t,y)

% SET GLOBAL VARIABLES
global A C n u
global s
global yb

% ADD THE SURFACE CONCENTRATION TO THE INTEGRATION ARRAY Y AND TEMPORARILY STORE AS Y1 ARRAY
y1 = [y' yb];

%--------------------------------------------------------
% THERE ARE TWO WAYS OF CALCULATING THE DERIVATIVE ARRAYS
% IF ID==1, USE THE COLLOCATION C MATRIX (THIS IS MUCH FASTER)
% IF ID==2, USE THE COLLOCATION A MATRIX
%--------------------------------------------------------
ID =2;

%-----------
% FOR ID = 1
%-----------
if ID==1
    for i=1:n;
        yprime(i) = sum(C(i,1:n+1).*y1);
    end
    
%-----------
% FOR ID = 2
%-----------
else
    for i=1:n
        sum1 	= 0;
        for j = 1:n+1
            sum2	= sum( A(j,1:n+1).*y1 );
            sum1	= sum1 + 2*u(j)^((1+s)/2)*A(i,j)*sum2;
        end
        yprime(i) 	= 2*u(i)^((1-s)/2)*sum1;
    end
    
end

yprime = yprime'; % Convert the derivative array into column format
