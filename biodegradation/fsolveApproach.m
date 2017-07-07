clear all; close all; clc;
global C T X c0 v
data = loadData('4NPfirst');
testNum=1;

pSize = length(data{1,testNum}.C);
indexes = round(linspace(1,pSize,pSize));
C = data{1,testNum}.C(indexes(2:end))';
T = data{1,testNum}.T(indexes(2:end))';
X = data{1,testNum}.X(indexes(2:end))';
c0 = data{1,testNum}.C(1);
v = data{1,testNum}.v;

[k, kAll] = findKmultistart('lsqnonlin', 1);
close all;
C=[c0;C];
T=[0;T];
plotCT(k);
plotRC(k);

function plotCT(k)
    global C T X c0 v
    subplot(2,1,1);
    ssv = mean(X);
    funK = @(T)C_runge(k, T, c0, ssv, v);
    hold on;
    C_interp = @(t)interp1(T, funK(T), t, 'linear');
    fplot(@(t)C_interp(t), [0  max(T)]);
    plot(T, funK(T), 'bo', 'LineWidth', 2);
    axis([0,inf,0,inf]);
    xlabel('t');
    ylabel('C');
    title(['K_{max}=' num2str(k(1)) ' K_s=' num2str(k(2)) ' K_I=' num2str(k(3))]);
    plot(T, C, 'rx', 'LineWidth', 2);
end

function plotRC(k)
    global C T X c0 v
    subplot(2,1,2);
    fplot(@(C)R_cal(k(2:3), C), [0  1.1*max(C)]);
    hold on;
    plot(C, R_cal(k(2:3), C), 'bo', 'LineWidth', 2);
    plot(k(2:3), R_cal(k(2:3), k(2:3)), 'bx', 'LineWidth', 2);
    axis([0,inf,0,inf]);
    title(num2str(R_cal(k(2:3), k(2:3))));
    xlabel('C');
    ylabel('R');
end

function [Ff, J] = implicitODE(k)
    global C T X c0 v
    
    Ff = 2*C.*k(3) + C.^2 + 2*k(2).*k(3).*log(C)-...
        -2*c0.*k(3) - c0.^2 - 2*k(2)*k(3)*log(c0) + 2.*(X./v)*k(1)*k(3).*T;
    
    J = [2*k(3).*(X./v).*T, ...
       2*k(3).*log(C) - 2*k(3).*log(c0),...
       2*C - 2*c0 + 2*k(2).*log(C) + 2*k(1).*(X./v).*T - 2*k(2)*log(c0)];
end

function [k, Kall] = findKmultistart(methodName, Nstart)
    global C
    startPoints = max(C)*randn(Nstart,3);
    Kall = zeros(Nstart,3);
    valBest = Inf;
    indBest = -1;
    for i = 1:Nstart
        x0 = startPoints(i,:);
    
        if strcmp(methodName, 'lsqnonlin')
            [k, cost] = solveLsqnonlin(x0);

        elseif strcmp(methodName, 'fsolve')
            [k, cost] = solveFsolve(x0);
        end
        
        Kall(i,:) = k;
        if norm(cost) <= valBest
            valBest = norm(cost);
            indBest = i;
        end
    end
    k = Kall(indBest,:);
end
function  [k, cost] = solveLsqnonlin(x0)
    global C T X c0 v
    options = optimoptions('lsqnonlin',...
        'Display','iter-detailed',...
        ...%'PlotFcn',@optimplotfirstorderopt,...
        'SpecifyObjectiveGradient',true,...
        'CheckGradients',true,'FiniteDifferenceStepSize',1e-4,...
        'Algorithm','trust-region-reflective',...%levenberg-marquardt
        'OptimalityTolerance',1e-14,'FunctionTolerance', 1e-14,...
        'MaxIterations',20000, 'MaxFunctionEvaluations', 60000);
    lb = [0,0,0];
    up = [max(C), max(C), max(C)];
    [k, cost] = lsqnonlin(@implicitODE, x0, lb, up, options);
end
function  [k, cost] = solveFsolve(x0)
    problem.options = optimoptions('fsolve',...
        'Display','iter-detailed',...
        'PlotFcn',@optimplotfirstorderopt,...
        'SpecifyObjectiveGradient',true,...
        'CheckGradients',true,'FiniteDifferenceStepSize',1e-5,...
        'Algorithm','levenberg-marquardt',...
        'OptimalityTolerance',1e-14,'FunctionTolerance', 1e-14,...
        'MaxIterations', 20000, 'MaxFunctionEvaluations', 60000);
    problem.objective = @implicitODE;
    problem.solver = 'fsolve';
    problem.x0 = x0;
    [k, cost] = fsolve(problem);
end

function y = C_runge(k, T, C0, ssv, v)
    %Runge Kutta
    options = odeset('RelTol',1e-4, 'AbsTol', 1e-4);
    [~,y] = ode23t(@(T, C)F(T, C), T, C0, options);
    y=y';

    function dCdt = F(t, c)
        dCdt = -k(1)*ssv*c/(v*(c+k(2)+c.^2/k(3)));
    end
end

function y = R_cal(k, C)
    y = C.*(1+2*sqrt(k(1)/k(2)))./(k(1)+C+C.^2/k(2));
end
