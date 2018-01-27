clear all; close all; clc;

load results_best_fmin_search.mat

global FSOLVE LSQNONLIN FMINCON FMINUNC FMINCON_INTERIOR LSQNONLIN_LM FMINUNC_QN FMINCON_SQP FMINCON_ACTIVESET
global TYPE_R TYPE_C TYPE_RUNGE
TYPE_C = 1; TYPE_RUNGE = 2; TYPE_R = 3;
FSOLVE = 1; LSQNONLIN = 2; LSQNONLIN_LM = 3; FMINUNC = 4; FMINUNC_QN = 5; FMINCON = 6; FMINCON_INTERIOR = 7; FMINCON_SQP = 8; FMINCON_ACTIVESET = 9;
ALGOSET = [FSOLVE LSQNONLIN LSQNONLIN_LM FMINUNC FMINUNC_QN FMINCON FMINCON_INTERIOR FMINCON_SQP FMINCON_ACTIVESET];

global isRegul lamdaReg
isRegul = true;
lamdaReg = 1e-10;
data = loadData();
testsSize = length(data);
lamdaSet = [0.01 0.03 0.1 0.3 0.5 1 5 10 30 100];
ALGORITHM = FMINCON_INTERIOR;
FUNCTION_TYPE = TYPE_C;
AUGMENTATION = true;
WITH_INEQ = false;

tables = cell(1, testsSize);
for testNum = 1:testsSize
    disp(['dataSet #' num2str(testNum)])
    pSize = length(data{testNum}.C);
    indexes = round(linspace(1,pSize,pSize));
    X = data{testNum}.X(indexes);
    T = data{testNum}.T(indexes);
    C = data{testNum}.C(indexes);
    ssv = mean(X);
    C0 = C(1);
    
    lb(1:3) = 1e-10;
    ub(1:3) = max(C);
    if WITH_INEQ
        A = zeros(3,3);
        A(2,2:3)=[1 -1];
        b = zeros(3,1);
    else
        A = [];
        b = [];
    end
    
    %kFirst = findMinimum(C, T, X, [], lb, ub, [], [], TYPE_RUNGE, LSQNONLIN_LM);
    
    len = length(fminBest{testNum}.Properties.RowNames);
    dataR = zeros(len,4);
    for i = 1:len
        kFirst = table2array(fminBest{testNum}(i, 1:3));
        data{testNum}.kFirst=kFirst;
        Cfit = dCdT(kFirst, T, C0, ssv, 1e-5);
        if length(C) == length(Cfit)
            checkSize = (size(C) ~= size(Cfit));
            if checkSize(1) || checkSize(2)
                Cfit = Cfit';
            end
            [r2,rsme]=rsquare(C, Cfit);
            dataR(i, 3) = r2;
        end
        if(AUGMENTATION == true)
            
            Taug = linspace(T(1), T(end), 101)';
            Caug = dCdT(kFirst, Taug, C0, ssv, 1e-4);
            %Caug = CfunFit(Taug, T, C);
            
            Raug = -diff(Caug(1:2:end))./diff(Taug(1:2:end));
            
            Raug = Raug ./ max(Raug);
            Taug = Taug(2:2:end);
            Caug = Caug(2:2:end);
            %A(2:3,2:3) b(2:3)
            kSecond = findMinimum(Caug, Taug, X, Raug, lb(1:2), ub(1:2), [] , [] , TYPE_R, FMINCON);
            dataR(i, 1:2) = kSecond;
            R = RfunFit(C(2:end), Caug, Raug);
        else
            R = -diff(C)./diff(T);
            CC = diff(C)/2 + C(1:end-1);
            %R = R ./ max(R);
            kSecond = findMinimum(C(2:end), T, X, R, lb(1:2), ub(1:2), A, b, TYPE_R, FMINCON);
        end
%         data{testNum}.R = R;
%         
% 
%         disp(['SecondApprox(k_s, k_I)=' num2str(kSecond)]);
%         
%         
%         
%         %cFunK = @(k)solveMassBalance23(k, T, C0, ssv, 1e-3, kSecond(1), kSecond(2));
%         kMax = kFirst(1);
%         %kMax = fitParameters(cFunK, C, lb(1), ub(1), [], [], 30);
%         data{testNum}.kSecond = [kMax kSecond];
%         disp(['SecondApprox(k_max)=' num2str(kMax)]);
%         %plotData(data);
    end
    
    dataR(:,4) = sqrt(dataR(:,1).*dataR(:,2));
    tt = array2table(dataR, 'VariableNames',{'K_s2','K_I2', 'r22', 'c_zv2'});
    tables{testNum}=[fminBest{testNum} tt];
end
save lastRun tables

%%
function k = findMinimum(C, T, X, R, lb, ub, A, b, FUNCTION_TYPE, ALGORITHM)
global FSOLVE LSQNONLIN FMINCON FMINUNC FMINCON_INTERIOR LSQNONLIN_LM FMINUNC_QN FMINCON_SQP FMINCON_ACTIVESET
global TYPE_R TYPE_C TYPE_RUNGE

tol = 1e-8;
iterations = 1000;
time = 2*60;
hasFunDeriv = true;
startPointIteration = round(iterations/3);

if FUNCTION_TYPE == TYPE_R
    funGoal = @(k)Rcal(k, C, R);
    funGoalSum = @(k)RcalSum(k, C, R);
elseif FUNCTION_TYPE == TYPE_C
    funGoal = @(k)Fsystem(k, C, T, X);
    funGoalSum = @(k)FsystemSum(k, C, T, X);
elseif FUNCTION_TYPE == TYPE_RUNGE
    ssv = mean(X);
    funGoal =  @(k)dCdT(k, T, C(1), ssv, sqrt(tol)) - C;
    funGoalSum =  @(k)sum(funGoal(k).^2);
    startPointIteration = round(startPointIteration/10);
    hasFunDeriv = false;
end
x0 = getStartPoint(funGoalSum, lb, ub, A, b, startPointIteration);

if ALGORITHM == FSOLVE
    optionsName = 'fsolve';
    options = getGeneralOptions(optionsName, tol, iterations);
    options.SpecifyObjectiveGradient = hasFunDeriv;
    options.Algorithm = 'trust-region-dogleg';
    k = fsolve(funGoal, x0, options);
elseif ALGORITHM == LSQNONLIN || ALGORITHM == LSQNONLIN_LM
    optionsName = 'lsqnonlin';
    options = getGeneralOptions(optionsName, tol, iterations);
    options.SpecifyObjectiveGradient = hasFunDeriv;
    if ALGORITHM == LSQNONLIN_LM
        options.Algorithm = 'levenberg-marquardt';
    else
        options.Algorithm = 'trust-region-reflective';
        if ~hasFunDeriv
            k = [0 0 0];
            return;
        end
    end
    k = multistart(optionsName, funGoal, x0, lb, ub, options, round(iterations / 10));
elseif ALGORITHM == FMINUNC || ALGORITHM == FMINUNC_QN
    optionsName = 'fminunc';
    options = getGeneralOptions(optionsName, tol, iterations);
    options.SpecifyObjectiveGradient = hasFunDeriv;
    if ALGORITHM == FMINUNC_QN
        options.Algorithm = 'quasi-newton';
    else
        options.Algorithm = 'trust-region';
    end
    if hasFunDeriv
        options.HessianFcn = 'objective';
    end
    %k = fminunc(funGoalSum, x0, options);
    k = multistart(optionsName, funGoalSum, x0, lb, ub, options, round(iterations / 10));
elseif ALGORITHM == FMINCON || ALGORITHM == FMINCON_INTERIOR || ...
        ALGORITHM == FMINCON_SQP || ALGORITHM == FMINCON_ACTIVESET
    optionsName = 'fmincon';
    options = getGeneralOptions(optionsName, tol, iterations);
    options.SpecifyObjectiveGradient = hasFunDeriv;
    if hasFunDeriv
        options.HessianFcn = @(x,lambda)objectiveToMultiHessian(x,lambda, funGoalSum);
    end
    if ALGORITHM == FMINCON_INTERIOR
        options.Algorithm = 'interior-point';
    elseif ALGORITHM == FMINCON_SQP
        options.Algorithm = 'sqp';
    elseif ALGORITHM == FMINCON_ACTIVESET
        options.Algorithm = 'active-set';
    else
        options.Algorithm = 'trust-region-reflective';
        if hasFunDeriv
            options.HessianFcn = 'objective';
        else
            k = [0 0 0];
            return;
        end
    end
    k = multiStartSolver(funGoalSum, x0, lb, ub, A, b, options, time);
end
end
%%
function k = multistart(optionsName, funGoal, x0, lb, ub, options, iterations)
Cproblem = createOptimProblem(optionsName, 'x0', x0,'objective', funGoal, 'options', options, 'lb', lb, 'ub', ub);
ms = MultiStart;
if(length(iterations) > 1)
    pts = CustomStartPointSet(iterations);
    %iterations = list(tpoints);
else
    pts = iterations;
end
k = run(ms, Cproblem, pts);
end
%% interpolation function for C(t) and R(c)
function Caug = CfunFit(t, T, C)
%C_interp = @(t)interp1(T, funK(T), t, 'linear');
[xData, yData] = prepareCurveData(T, C);
%'concavedown', [0, max(C)*0.3],
fitresultC = slmengine(xData, yData, 'plot', 'off', 'decreasing', 'on','rightminvalue', 0);
Caug=slmeval(t, fitresultC);
end
%%
function R = RfunFit(c, Caug, Raug)
[xData, yData] = prepareCurveData(Caug, Raug);
ft = fittype('smoothingspline');
opts = fitoptions('Method', 'SmoothingSpline', 'SmoothingParam', 0.9999);
[Rinter, gof] = fit(xData, yData, ft, opts);
R = Rinter(c);
end
%% evaluation
function [r2,rmse] = evaluate(Yexp, Yfit)
[r2, rmse] = rsquare(Yexp, Yfit);
end
%% Approach I. Get a nonlinear system or least square function after implicite solution of ODE in each point
%  Use: fslove 'trust-region-dogled', lsqnonlin, fmincon, fmincon
function [y, grad, H] = FsystemSum(k, C, T, X)
global lamdaReg
[fun, J] = Fsystem(k, C, T, X);
fun = fun + (lamdaReg/2)*sum(k.^2);
J = J + lamdaReg*k;
y=sum(fun.^2);
if nargout > 1
    gradi = 2*fun.*J;
    grad=sum(gradi);
    
    if nargout > 2
        N = length(fun);
        C0 = C(1);
        X = X(2:end);
        C = C(2:end);
        T = T(2:end);
        
        HI = zeros(3,3,N);
        HI(1,3, :) = 2*X.*T;
        HI(3,1, :) = HI(1,3, :);
        HI(2,3, :) = 2*log(C) - 2*log(C0);
        HI(3,2, :) = HI(2,3, :);
        
        Hsum =  zeros(3,3);
        for i = 1:N
            Hsum = Hsum + fun(i).*HI(:,:, i);
        end
        H = 2*((gradi')*gradi + Hsum) + lamdaReg;
    end
end
end
%%
function [F, J] = Fsystem(k, C, T, X)
C0 = C(1);
X = X(2:end);
C = C(2:end);
T = T(2:end);

F = 2*X.*k(1).*k(3).*T + 2*k(2).*k(3).*log(C) + 2*C.*k(3) + C.^2 -...
    -2*C0.*k(3) - C0.^2 - 2*k(2)*k(3)*log(C0);
if nargout > 1
    J = [2*k(3).*X.*T, ...
        2*k(3).*log(C) - 2*k(3).*log(C0),...
        2*k(1).*X.*T  + 2*k(2).*log(C) + 2*C - 2*C0  - 2*k(2)*log(C0)];
end
end
%%
function H = objectiveToMultiHessian(k, lambda, fun)
[~, ~, H] = fun(k');

% k2 = k(2); k3 = k(3);
% HC = zeros(3,3);
% HC(2:3,2:3) = [-k3^2./(4*(k2.*k3).^(3/2)), 1./(2*(k2.*k3).^(1/2)) - (k2.*k3)./(4*(k2.*k3).^(3/2));
%     1./(2*(k2.*k3).^(1/2)) - (k2.*k3)/(4*(k2.*k3).^(3/2)), -k2.^2/(4.*(k2.*k3).^(3/2))];
%HEND = H + lambda.ineqnonlin(1)*HC - lambda.ineqnonlin(2)*HC;

end
%% Tomei approach using R
function y = R_cal(k, C)
% syms c k1 k2
% rsym = c*(1+2*sqrt(k1/k2))/(k1+c+c^2/k2);
% dk1 = diff(rsym, k1);
% dk2 = diff(rsym, k2);
% dk1k1 = simplify(diff(dk1, k1))
% dk1k2 = simplify(diff(dk1, k2))
% dk2k2 = simplify(diff(dk2, k2))
y = C.*(1+2*sqrt(k(1)/k(2)))./(k(1)+C+C.^2/k(2));
end
%%
function [y, J] = Rcal(k, C, R)
y = R_cal(k, C) - R;
if nargout > 1
    J(:,1)=C./(k(2).*sqrt(k(1)./k(2))*(C.^2./k(2) + C + k(1))) - (C.*(2*sqrt(k(1)./k(2)) + 1))./(C + k(1) + C.^2/k(2)).^2;
    J(:,2)=(C.*(C.^2.*k(1) - k(1).^2.*k(2) + C.^2.*k(2).*sqrt(k(1)./k(2)) - C.*k(1).*k(2)))./(k(2).*sqrt(k(1)./k(2))*(C.^2 + k(2)*C + k(1)*k(2)).^2);
end
end
%%
function [y, grad, H] = RcalSum(k, C, R)
global lamdaReg
if nargout > 1
    [fun, J] = Rcal(k, C, R);
    J = J + lamdaReg*k;
else
    fun = Rcal(k, C, R);
end
fun = fun + (lamdaReg/2)*sum(k.^2);
y=sum(fun.^2);
if nargout > 1
    gradi = 2*fun.*J;
    grad=sum(gradi);
    if nargout > 2
        if nargout > 2
            N = length(fun);
            
            HI = zeros(2,2,N);
            HI(1,1, :) = -(C.*(2*C.^3*k(1) - 4*k(1).^4*(k(1)/k(1)).^(3/2) + C.^4 + C.^2*k(1).^2 - 3*k(1).^2*k(1).^2 + 6*C.*k(1)*k(1).^2 + 6*C.^2*k(1)*k(1)))./(2*k(1)*(k(1)/k(1)).^(3/2)*(C.^2 + k(1)*C + k(1)*k(1)).^3);
            HI(1,2, :) = -(C.^3*k(1).^2 - C.^5 - C.*k(1).^2*k(1).^2 + 4*C.^3*k(1).^2*(k(1)/k(1)).^(1/2) + 6*C.^3*k(1)*k(1))./(2*k(1)*(k(1)/k(1)).^(1/2)*(C.^2 + k(1)*C + k(1)*k(1)).^3);
            HI(2,1, :) = HI(1,2, :);
            HI(2,2, :) = -(C.*(C.^4*k(1).^2 - 3*k(1).^4*k(1).^2 - 6*C.*k(1).^3*k(1).^2 + 6*C.^2*k(1).^3*k(1) + 6*C.^3*k(1).^2*k(1) - 3*C.^2*k(1).^2*k(1).^2 + 4*C.^3*k(1).^3*(k(1)/k(1)).^(3/2) + 4*C.^2*k(1)*k(1).^3*(k(1)/k(1)).^(3/2)))./(2*k(1).^3*(k(1)/k(1)).^(3/2)*(C.^2 + k(1)*C + k(1)*k(1)).^3);
            
            Hsum =  zeros(2,2);
            for i = 1:N
                Hsum = Hsum + fun(i).*HI(:,:, i);
            end
            H = 2*((gradi')*gradi + Hsum)+lamdaReg;
        end
    end
end
end
%%
function [x0, population] = getStartPoint(funGoal, lb, ub, A, b, generations)
nvars = length(lb);
disp('find start point x0 using genetic random search...');
optionsGA = gaoptimset('PopInitRange',  [min(lb); max(ub)], 'Display', 'iter',...
    'PopulationSize', 100, 'Generations', generations);
[x0,~,~,~,population,~] = ga(funGoal, nvars, A, b, [], [], lb, ub, [], [], optionsGA);
%x0 = (ub - lb)./2;
end
%%
function [c,ceq] = nonlcon(k)
c = [sqrt(k(2)*k(3));-sqrt(k(2)*k(3))];
ceq = [k(3); -k(2)];
end
%%
function k = multiStartSolver(funGoal, x0, lb, ub, A, b, options, time)
disp('MultiSearch running...');
Cproblem = createOptimProblem('fmincon', 'x0', x0, 'objective', funGoal,...
    'Aineq', A, 'bineq', b, 'lb', lb, 'ub', ub, 'options', options);
%Cproblem.nonlcon = @(k)nonlcon(k);
k = run(GlobalSearch('StartPointsToRun','bounds-ineqs',...
    'MaxWaitCycle', ceil(options.MaxIterations/50),'NumStageOnePoints', ceil(options.MaxIterations/5), ...
    'NumTrialPoints', options.MaxIterations, 'MaxTime', time), Cproblem);
end
%%
function c = dCdT(k, T, C0, ssv, tol)
options = odeset('Nonnegative', 1, 'RelTol', tol, 'AbsTol', tol*0.1);
%ode45, ode15s, ode23t
%display(k)
% if k(2) < 1e-4
%      k(2)= 1e-4;
% end
[~,c] = ode15s(@(t,c)-k(1)*ssv.*c/(k(2)+c+c.^2./k(3)), T, C0, options);
if length(c) ~= length(T)
    c = C_runge(k, T, C0, ssv);
end
end

%%
function y = C_runge(k, T, C0, ssv)
%Euler method
size=length(T);
y=zeros(1,size);
for i=1:size
    if(i == 1)
        y(1)=C0;
    else
        y(i) = y(i-1) + F(0, y(i-1)) * (T(i)-T(i-1));
    end
end

    function dCdt = F(t, c)
        dCdt = -k(1)*ssv*c/((c+k(2)+c.^2/k(3)));
    end
end
%%
function c = solveMassBalance23(k, T, C0, ssv, tol, k2, k3)
options = odeset('Nonnegative', 1, 'RelTol', tol, 'AbsTol', tol*0.1);
[~,c] = ode23t(@(t,c)-k(1)*ssv.*c/(k2+c+c.^2./k3), T, C0, options);
end
%%
function options = getGeneralOptions(optionsName, tol, iterations)
options = optimoptions(optionsName,...
    'Display', 'iter-detailed',...
    'OptimalityTolerance', tol, 'FunctionTolerance', tol,...
    'MaxIterations', iterations, 'MaxFunctionEvaluations', 6*iterations);
%options.PlotFcn = @optimplotfirstorderopt;
%options.FiniteDifferenceStepSize = 1e-12;
%options.FiniteDifferenceType = 'central';
%options.CheckGradients = true;
%options.UseParallel = true;
end
%%
function k = fitParameters(funK, funExp, lb, ub, A, b)
funGoal =  @(k)sum((funK(k)-funExp).^2);
x0 = getStartPoint(funGoal, lb, ub, A, b, 100);
k = multiStartSolver(funGoal, x0, lb, ub, A, b, 1e-7, 1000, 2*60);
end
%% plot data
function plotData(data)
figure();
testsSize = 1;%length(data);
for testNum = 1:testsSize
    for frameInd = 1:3
        titleD = data{testNum}.title;
        X = data{testNum}.X;
        T = data{testNum}.T;
        C = data{testNum}.C;
        R = data{testNum}.R;
        ssv = mean(X);
        C0 = C(1);
        
        if frameInd == 1
            k = data{testNum}.kFirst;
        else
            k = data{testNum}.kSecond;
        end
        
        subplot(testsSize, 3, (testNum - 1) * 3 + frameInd);
        hold on;
        if frameInd == 2
            Cmodel = linspace(0, max(C), 100);
            Rmodel = R_cal(k(2:3), Cmodel);
            
            plot(Cmodel, Rmodel, 'b-');
            plot(C(2:end), R, 'rx', 'LineWidth', 2);
            
            Rfit = R_cal(k(2:3), C(2:end));
            r2=rsquare(R, Rfit);
            title({['[' titleD '] K_s=' num2str(k(2)) ' K_I=' num2str(k(3))], ['r2=' num2str(r2)]});
            xlabel('C');
            ylabel('R');
            axis([0,inf,0,inf]);
        elseif frameInd == 1 || frameInd == 3
            Tmodel=linspace(0, max(T), 100);
            Cmodel = dCdT(k, Tmodel, C0, ssv, 1e-4);
            
            plot(Tmodel, Cmodel, 'b-');
            plot(T, C, 'rx', 'LineWidth', 2);
            
            Cfit = dCdT(k, T, C0, ssv, 1e-4);
            r2=rsquare(C, Cfit);
            title({['[' titleD '] K_{max}=' num2str(k(1)) ' K_s=' num2str(k(2)) ' K_I=' num2str(k(3))], ['r2=' num2str(r2)]});
            xlabel('t');
            ylabel('C');
            axis([0, inf, 0, inf]);
        end
    end
end
end
%%
function findForAllAlgorithms()
FUNSET = [TYPE_C, TYPE_RUNGE];
kFirst = zeros(length(FUNSET), 3, length(ALGOSET));
kFirstGOF  = zeros(length(FUNSET), 2, length(ALGOSET));
for fun=FUNSET
    %algo = FMINCON;
    for algo=ALGOSET(FMINCON_INTERIOR)
        k = findMinimum(C, T, X, [], lb, ub, A, b, fun, algo);
        Cfit = dCdT(k, T, C0, ssv, 1e-4);
        
        if length(C) == length(Cfit)
            [r2,rsme]=rsquare(C, Cfit);
            kFirstGOF(fun, :, algo) = [r2,rsme];
        end
        kFirst(fun, :, algo) = k;
    end
end
testRes{testNum}.kFirst = kFirst;
testRes{testNum}.kFirstGOF = kFirstGOF;
end











