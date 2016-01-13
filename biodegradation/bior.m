%C_NP = c;
%k_max = k1;
%K_s = k2;
%K_I = k3;
%================================
%dc/dt = (k1*X/V)*c/(c+k2+c^2/k3);
%Solution t:
%t = V/(k1*X*k3)*(k2*k3*log(c)+k3*c+c^2/2) + const
%x(1) = V/(k1*X*k3)*(k2*k3*log(x(2))+k3*x(2)+x(2)^2/2) + x(3)
%Solution c:
%c_cal(t) = Integral(0, ti, k3*k1*X(t)*c(t)/(V*k3*(c(t)+k2+c(t)^2)))dt

clear;
close all;
data = loadData();
m=length(data);

nvars = 2;
R=[];
C=[];
X=[];
for testNum = 1:m
    T = data{testNum}.T;
    Cr = data{testNum}.C;
    R = [R -diff(Cr)./diff(T)];
    
    C=[C Cr(1:end-1)];
    Xr = data{testNum}.X;
    X=[X Xr(1:end-1)];
end
v = data{1}.v; %l, litres
n = length(T);
ssv = mean(X);
R_exp = R ./ max(R);
 
SORT = [R;C;X];
SORT = sort(SORT, 2);
R=SORT(1,:);
C=SORT(2,:);
X=SORT(3,:);

lb(1:2) = 1e-10;
ub(1:2) = max(C);
x0 = (ub-lb)./2;
A = [1 -1]; b = 0;
R_cal = @(k, C)C.*(1+2*sqrt(k(1)/k(2)))./(k(1)+C+C.^2/k(2));
main_fun_R = @(k)sum((R_cal(k, C) - R_exp).^2);
optionsGA = gaoptimset('PopInitRange',  [min(lb); max(ub)], 'Display', 'iter');
optionsGA.PopulationSize = 200;
optionsGA.Generations = 20;
x0 = ga(main_fun_R, nvars,A,b,[],[],lb,ub,[],[],optionsGA);

[x0 errfun] = ga(main_fun_R, nvars,[],[],[],[],lb,ub,[],[],optionsGA)
Rproblem = createOptimProblem('fmincon', 'x0', x0(1:nvars),...
    'objective', main_fun_R,'Aineq', A, 'bineq', b, 'lb', lb(1:nvars), 'ub', ub(1:nvars));
[RK,errormulti] = run(GlobalSearch,Rproblem)

data{testNum}.RK = RK;

close all;
fplot(@(C)R_cal(RK, C), [0  1.1*max(C)]);
hold on;
plot(C, R_cal(RK, C), 'bo', 'LineWidth', 2);
axis([0,inf,0,inf]);
xlabel('C');
ylabel('R');
title(['[' data{testNum}.title '] K_{s}=' num2str(RK(1)) ' K_I=' num2str(RK(2)) ]);
plot(C, R_exp, 'rx', 'LineWidth', 2);



