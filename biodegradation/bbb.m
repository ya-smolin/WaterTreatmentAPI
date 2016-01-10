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
%global v;
v = 2.25; %l, litres
dimensionless = false;
for testNum = 1:m
    if testNum ~= 1 continue; end;
    %global X T C n;
    
    X = data{testNum}.X;
    x=mean(X);
    T = data{testNum}.T;
    C = data{testNum}.C;
    n = length(T);
    
    %t = V/(k1*X*k3)*(k2*k3*log(c)+k3*c+c^2/2) + const
    
    %F = @(k) sum(abs(-T - v./(k(1).*x.*k(3)) * (k(2)*k(3)*log(C)+k(3)*C+C.^2/2)+k(4)));
    b=1;
    e=length(C);
    c1 = C(b:e-1);
    c2 = C(b+1:e);
    t1 = T(b:e-1);
    t2=T(b+1:e);
    x1 = X(b:e-1);
    x2 = X(b+1:e);
    f_int = @(k, c)(v./x).*(k(2)*k(3).*log(c)+k(3).*c+c.^2/2)./(k(1)*k(3));
    Fsystem = @(k) t2 - t1 + (f_int(k, c1) - f_int(k, c2));
    Ffun = @(k)sum((Fsystem(k)).^2);
    
    
    nvars = 3;
    lb(1:3) = 1e-10;
    ub(1:3) = 5*max(C); %1e10;
    generSize = 150;
    populSize = 100;
    optionsGA=gaoptimset();
    optionsGA = gaoptimset(optionsGA,'PopInitRange',  [min(lb); max(ub)],'PopulationSize', populSize,...
        'Generations', generSize,'Display', 'iter','TolFun',1e-12);
%     [x0 funV] = ga(Ffun, nvars,[],[],[],[],lb,ub,[],[],optionsGA)
%     disp(x0);
    

    p0 = 5*ones(1,3);
    problem = createOptimProblem('lsqnonlin','x0',p0,'objective',Fsystem,...
        'lb',lb,'ub',ub);
    ms = MultiStart('PlotFcns',@gsplotbestf);
    [xmulti,errormulti] = run(ms,problem,5000)

  
%     optionsFit = optimoptions('lsqnonlin');
%     optionsFit = optimoptions(optionsFit,'Display', 'iter', 'MaxFunEvals', 3000, 'TolX',  1e-16, 'TolFun', 1e-16);
%     [x Xfun] = lsqnonlin(Fsystem,x0,lb,ub, optionsFit)
    %K = fsolve(Fsystem, x0, options);
    
 
%     data{testNum}.RK = RK;
%       
%     subplot(m,2,2*testNum-1);
%     hold on;
%     fplot(@(C)R_cal(RK, C), [0  1.1*max(C)]);
%     plot(C, R_cal(RK, C), 'bo', 'LineWidth', 2);
%     axis([0,inf,0,inf]);
%     xlabel('C');
%     ylabel('R');
%     title(['[' data{testNum}.title '] K_{s}=' num2str(RK(1)) ' K_I=' num2str(RK(2)) ]);
%     
%     plot(C, R_exp, 'rx', 'LineWidth', 2);
%   
%     subplot(m,2,2*testNum);
%     hold on;
%     C_interp = @(t)interp1(T,C_cal(CK, T),t, 'spline');
%     fplot(@(t)C_interp(t), [0  max(T)]);
%     plot(T, C_cal(CK, T), 'bo', 'LineWidth', 2);
%     axis([0,inf,0,inf]);
%     xlabel('t');
%     ylabel('C');
%     title(['[' data{testNum}.title '] K_{max}=' num2str(CK(1)) ' K_s=' num2str(CK(2)) ' K_I=' num2str(CK(3)) ]);
%     
%     plot(T, C, 'rx', 'LineWidth', 2);
    
end