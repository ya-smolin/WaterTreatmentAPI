c = [91.5
78
74.5
56
42.5
10
1
0.3];

t = [0
26.4
48
69.6
81.6
112.8
124.8
144];

X = [205
203
186
188
189
184
210
187
];
V = 2.2; %litres

%C_NP = c;
%k_max = k1;
%K_s = k2;
%K_I = k3;
%================================
%dc/dt = (k1*X/V)*c/(c+k2+c^2/k3);
%Solution t:
%t = V/(k1*X*k3)*(k2*k3*log(c)+k3*c+c^2/2) + const
%Solution c:
%c_cal(t) = Integral(0, ti, k3*k1*X(t)*c(t)/(V*k3*(c(t)+k2+c(t)^2)))dt
%
iif = @(varargin) varargin{2*find([varargin{1:2:end}], 1, 'first')}();

F = @(k, i) k(3)*k(1)*X(i)*c(i)/(V*(k(3)*c(i)+k(3)*k(2)+c(i).^2));

c_calc_rec = @(k, i, frec) iif(i==1,c(1),...
    true, @()frec(k, i-1, frec) + (F(k, i)+F(k, i-1))/2*(t(i)-t(i-1)));
c_calc = @(k, i)c_calc_rec(k, i, c_calc_rec);

%F(1,2,3,5)
main_fun_rec = @(k, i, frec) iif(i==0, 0,...
    true, @()frec(k, i-1,frec) + (c_calc(k, i) - c(i)).^2 );
Imax = length(t);
main_fun = @(k)main_fun_rec(k, Imax, main_fun_rec);


Cmax = max(c);
fun = main_fun;
nvars = 3;
lb = [0 0 0]';
ub =  [Cmax Cmax Cmax];
generSize = 5;
populSize = 50;

%% Start with the default options
options = gaoptimset;
%% Modify options setting
%options = gaoptimset(options,'FitnessLimit', 0.8);
%options = gaoptimset(options,'InitialPopulation', initialPopulation);
options = gaoptimset(options,'PopInitRange',  [0; Cmax]);
options = gaoptimset(options,'PopulationSize', populSize);
options = gaoptimset(options,'Generations', generSize);
%options = gaoptimset(options,'StallGenLimit', 20);
%options = gaoptimset(options,'TolFun', 0.0001);
options = gaoptimset(options,'Display', 'iter');
%options = gaoptimset(options,'EliteCount', 3);
%options = gaoptimset(options,'PlotFcns', { @gaplotbestf });
k0 = ga(fun, nvars,[],[],[],[],lb,ub,[],[],options);
disp(k0);


%% Start with the default options
options = optimoptions('fmincon');
%% Modify options setting
options = optimoptions(options,'Display', 'iter');
[k,fval,exitflag,output,lambda,grad,hessian] = ...
fmincon(main_fun,k0,[],[],[],[],lb,ub,[],options);
disp(k);
for i = 1:Imax
    c_calcAr(i) = c_calc(k, i);
end
plot(t, c_calcAr, 'b+');
hold on;
plot(t, c, 'g*');
