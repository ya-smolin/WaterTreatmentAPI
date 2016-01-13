clear;
close all;
data = loadData();

m = length(data);
% m=6;

nvars = 3;
lb(1:nvars) = 1e-10;
ub(1:nvars) = 1000;
x0 = (ub - lb)./2;
optionsGA = gaoptimset('PopInitRange',  [min(lb); max(ub)], 'Display', 'iter');
optionsGA.PopulationSize = 200;
optionsGA.Generations = 60;



%     C_fun = @(k, T) C_runge(k, T, C0, ssv, v);
main_fun_C = @(k)C_all(k, data);
A= [0 1 -1]; b = 0;
[x0 errfun] = ga(main_fun_C, nvars,A,b,[],[],lb,ub,[],[],optionsGA)
Cproblem = createOptimProblem('fmincon', 'x0', x0,'objective', main_fun_C,...
    'Aineq', A, 'bineq', b, 'lb',lb,'ub',ub);
[kC,errormulti] = run(GlobalSearch,Cproblem)


for testNum = 1:m
    T = data{testNum}.T;
    %%%%%%%%%%%%PLOT%%%%%%%%%%%%%%%%%%%%%%
    subplot(m, 1, testNum);
    hold on;
    T = data{testNum}.T;
    C = data{testNum}.C;
    v = data{testNum}.v;
    X = data{testNum}.v;
    n = length(T);
    ssv = mean(X);
    C_fun = @(k, T) C_runge(k, T, C(1), ssv, v);
    C_interp = @(t)interp1(T, C_fun(kC, T), t, 'linear');
    fplot(@(t)C_interp(t), [0  max(T)]);
    plot(T, C_fun(kC, T), 'bo', 'LineWidth', 2);
    axis([0,inf,0,inf]);
    xlabel('t');
    ylabel('C');
    title(['[' data{testNum}.title '] K_{max}=' num2str(kC(1)) ' K_s=' num2str(kC(2)) ' K_I=' num2str(kC(3)) ]);
    plot(T, C, 'rx', 'LineWidth', 2);
    
end