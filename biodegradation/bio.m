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
%t = V/(k1*X*k3)*(k2*k3*log(c)+k3*c+c^2/2) + const

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
global v;
v = 2.25; %l, litres
dimensionless = false;
for testNum = 1:m
%     if testNum ~= 1 continue; end;
    global X T C n;
    
    X = data{testNum}.X;
    %X(:) = mean(X);
    T = data{testNum}.T;
    C = data{testNum}.C;
    
    R = -diff(C)./diff(T);
    R(end+1) = R(end);
    R_exp = R ./max(R);
    n = length(T);
    
    %k1>0 k2>0
     R_cal = @(k, C)C.*(1+2*sqrt(k(1)/k(2)))./(k(1)+C+C.^2/k(2));
%     R_cal = @(k, C)C.*(1+2*sqrt(k(1)/k(2)))./(k(1)+C);
    
    if(dimensionless)
        X = X / max(C); %[0 .. 1]
        T = T / max(T); %[0 .. 1]
        C = C / max(C); %[0 .. 1]
        v = v / max(C); %litres
    end

    nvars = 3;
    lb(1:nvars) = 1e-10;
    ub(1:nvars) = 5*max(C);
    
    generSize = 20;
    populSize = 200;
    optionsGA=gaoptimset();
    optionsGA = gaoptimset(optionsGA,'PopInitRange',  [min(lb); max(ub)],'PopulationSize', populSize,...
        'Generations', generSize,'Display', 'iter');
    optionsFit = optimoptions('lsqcurvefit');
    %'Algorithm', 'trust-region-reflective'
    optionsFit = optimoptions(optionsFit,'Display', 'iter', 'MaxFunEvals', 3000, 'TolFun', 1e-8);
    
    nvars=2;
    main_fun_R = @(k)sum((R_cal(k, C)- R_exp).^2);
    Rk0 = ga(main_fun_R, nvars,[],[],[],[],lb(1:nvars),ub(1:nvars),[],[],optionsGA);
    RK = lsqcurvefit(R_cal,Rk0,C,R,lb(1:nvars),ub(1:nvars),optionsFit);
%     RK = fmincon(main_fun_R, Rk0,[],[],[],[],lb(1:nvars),ub(1:nvars));

    nvars=3;
    main_fun_C = @(k)sum((C_cal(k, T)- C).^2);
    Ck0 = ga(main_fun_C, nvars,[],[],[],[],lb(1:nvars),ub(1:nvars),[],[],optionsGA);
    CK = lsqcurvefit(@C_cal,Ck0,T,C,lb(1:nvars),ub(1:nvars),optionsFit);

    if(dimensionless)
        data{testNum}.CK = CK  * max(data{testNum}.C);
        X = X * max(data{testNum}.C); %[0 .. 1]
        T = T * max(data{testNum}.T); %[0 .. 1]
        C = C * max(data{testNum}.C); %[0 .. 1]
        v = v * max(data{testNum}.C); %litres
    else
        data{testNum}.CK = RK;
    end
    data{testNum}.RK = RK;
      
    subplot(m,2,2*testNum-1);
    hold on;
    fplot(@(C)R_cal(RK, C), [0  1.1*max(C)]);
    plot(C, R_cal(RK, C), 'bo', 'LineWidth', 2);
    axis([0,inf,0,inf]);
    xlabel('C');
    ylabel('R');
    title(['[' data{testNum}.title '] K_{s}=' num2str(RK(1)) ' K_I=' num2str(RK(2)) ]);
    
    plot(C, R_exp, 'rx', 'LineWidth', 2);
  
    subplot(m,2,2*testNum);
    hold on;
    C_interp = @(t)interp1(T,C_cal(CK, T),t, 'spline');
    fplot(@(t)C_interp(t), [0  max(T)]);
    plot(T, C_cal(CK, T), 'bo', 'LineWidth', 2);
    axis([0,inf,0,inf]);
    xlabel('t');
    ylabel('C');
    title(['[' data{testNum}.title '] K_{max}=' num2str(CK(1)) ' K_s=' num2str(CK(2)) ' K_I=' num2str(CK(3)) ]);
    
    plot(T, C, 'rx', 'LineWidth', 2);
    
end