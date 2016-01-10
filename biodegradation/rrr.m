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
test;
data{1}.C = [91.5 78 74.5 56 42.5 10 1 0.3]; %C, mg/l
data{1}.T = [0 26.4 48 69.6 81.6 112.8 124.8 144]; %Ò, hours
data{1}.X = [205 203 186 188 189 184 210 187]; %ss, mg/l
data{1}.title='90';

data{2}.C = [80.1 69.5 65.2 55.5 22 0.2]; %C, mg/l
data{2}.T = [0 14.88 24 57.6 96 115.2]; %Ò, hours
data{2}.X = [196 195 199 213 248 248]; %ss, mg/l
data{2}.title='80';

data{3}.C = [54.5 54 50.8 39 30 9 0]; %C, mg/l
data{3}.T = [0 6.96 24 54.96 69.6 96 110.4]; %Ò, hours
data{3}.X = [186 183 197 198 218 227 227]; %ss, mg/l
data{3}.title='60';

data{4}.C = [36.3 35.5 31 25.1 13.8 4 0]; %C, mg/l
data{4}.T = [0 0.167 1 2 3 3.5 4]*24; %Ò, hours
data{4}.X = [207 204 225 230 240 246 253]; %ss, mg/l
data{4}.title='40';

data{5}.C = [23.8 21.7 20 11 3 0.1]; %C, mg/l
data{5}.T = [0 0.8 1 2.5 3.3 4]*24; %Ò, hours
data{5}.X = [234 237 249 264 263 262]; %ss, mg/l
data{5}.title='20';

data{6}.C = [24.9 20.5 18.8 14 7.7 0.5 0]; %C, mg/l
data{6}.T = [0 18 24 43 55 72 91]; %Ò, hours
data{6}.X = [197 174 183 189 190 202 220]; %ss, mg/l
data{6}.title='20';
m=length(data);
global v;
v = 2.25; %l, litres
dimensionless = false;
isCalculateR = false;


global X T C ssv n;

lb(1:2) = 1e-10;
ub(1:2) = max(C);
x0 = (ub-lb)./2;

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

n = length(T);
ssv = mean(X);

C = Y(1:2:end)' + sign(rand(1) - 0.5)*rand(1);
T = tt(1:2:end)';
v=2.1;
ssv = 825;
R = -diff(C)./diff(T);
C = C(1:end-1);
T = T(1:end-1);
% k = [4.65 17.1 31.4];


R_exp = R ./ max(R);
 
% A = [R;C;X];
% A = sort(A, 2);
% R=A(1,:);
% C=A(2,:);
% X=A(3,:);

R_cal = @(k, C)C.*(1+2*sqrt(k(1)/k(2)))./(k(1)+C+C.^2/k(2));
main_fun_R = @(k)sum((R_cal(k, C) - R_exp).^2);
optionsGA = gaoptimset('PopInitRange',  [min(lb); max(ub)], 'Display', 'iter');
optionsGA.PopulationSize = 200;
optionsGA.Generations = 20;
x0 = ga(main_fun_R, nvars,[],[],[],[],lb,ub,[],[],optionsGA);

[x0 errfun] = ga(main_fun_R, nvars,[],[],[],[],lb,ub,[],[],optionsGA)
Rproblem = createOptimProblem('fmincon', 'x0', x0(1:nvars),...
    'objective', main_fun_R,...
    'lb', lb(1:nvars), 'ub', ub(1:nvars));
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



