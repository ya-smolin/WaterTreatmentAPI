% The inhibitor constant Ki is the concentration of inhibitor which is required to decrease the maximal rate of the reaction to half of the uninhibited value, in the presence of a low substrate concentration.

% Therefore, the lower the Ki the lower the concentration of inhibitor needed to lower the rate.

%#########################ANONYMOUS FUNCTION RECURSION####################
% iif = @(varargin) varargin{2*find([varargin{1:2:end}], 1, 'first')}();
% F = @(k, i) k(3)*k(1)*X(i)*C(i)/(v*(k(3)*C(i)+k(3)*k(2)+C(i).^2));
% c_calc_rec = @(k, i, frec) iif(i==1,C(1),...
%     true, @()frec(k, i-1, frec) - (F(k, i)+F(k, i-1))/2*(T(i)-T(i-1)));
% c_calc = @(k, i)c_calc_rec(k, i, c_calc_rec);
% main_fun_rec = @(k, i, frec) iif(i==0, 0,...
%     true, @()frec(k, i-1,frec) + (c_calc(k, i) - C(i)).^2 );
% main_fun = @(k)main_fun_rec(k, n, main_fun_rec);
% main_fun([1 2 3])

% #####################ALL IN ONE#########################################
% X=[];
% T=[];
% C=[];
% for i = 1:length(data)
%     X = [X data{i}.X];
%     %     X(:) = mean(X);
%     T = [T data{i}.T];
%     C = [C data{i}.C];
% end
% A = [T;C;X];
% A = sort(A, 1);
% T=A(1,:);
% C=A(2,:);
% X=A(3,:);


generSize = 20;
populSize = 200;
optionsGA=gaoptimset();
optionsGA = gaoptimset(optionsGA,'PopInitRange',  [min(lb); max(ub)],'PopulationSize', populSize,...
    'Generations', generSize,'Display', 'off');

 nvars=2;
 main_fun_R = @(k)sum((R_cal(k, C) - R_exp).^2);
 Rk0 = ga(main_fun_R, nvars,[],[],[],[],lb,ub,[],[],optionsGA);
 
 %     if(dimensionless)
%         X = X / max(C); %[0 .. 1]
%         T = T / max(T); %[0 .. 1]
%         C = C / max(C); %[0 .. 1]
%         v = v / max(C); %litres
%     end

%     if(dimensionless)
%         CK = CK  * max(data{testNum}.C);
%         X = X * max(data{testNum}.C); %[0 .. 1]
%         T = T * max(data{testNum}.T); %[0 .. 1]
%         C = C * max(data{testNum}.C); %[0 .. 1]
%         v = v * max(data{testNum}.C); %litres
%     end

%     [x0 errfun] = patternsearch(main_fun_C, x0, [],[],[],[],lb,ub, [], psoptimset('Display', 'iter', 'MaxIter', 1000))
%     [x0 errfun] =  particleswarm(main_fun_C,nvars,lb,ub,optimoptions(@particleswarm, 'Display', 'iter', 'StallIterLimit', 1000))
%     Cproblem = createOptimProblem('lsqcurvefit', 'x0', x0(1:nvars),'objective', @C_runge,...
%         'lb',lb(1:nvars),'ub',ub(1:nvars),'xdata', T, 'ydata', C);
%     [CK,errormulti] = run(ms,Cproblem,200) %disp


%     size=length(T);
%     y=zeros(1,size);
%     for i=1:size
%         if(i == 1)
%             y(1)=C(1);
%         else
%             %Real runge kutta
%             y(i) = y(i-1) + F(y(i-1), k, ssv, v) * (T(i)-T(i-1));
% 
%         end
%     end


  frequency = 455;
    T_ext = linspace(0, max(T), frequency);
    C_ext = C_runge(kFirst, T_ext, C0, ssv, v, false);
    R = -diff(C_ext(1:2:end))./diff(T_ext(1:2:end));
    T = T_ext(2:2:end);
    C = C_ext(2:2:end);
    
    
    function  [k, cost] = fitParametersFsolve(funGoal, lb, ub, A, b)
    nvars = length(lb);
    %x0 = (ub - lb)./2;
    %funGoal =  @(k)sum(funK(k).^2);
    
    disp('find the x_0 for the best approximation to the MultiSearch...');
    optionsGA = gaoptimset('PopInitRange',  [min(lb); max(ub)], 'Display', 'iter',...
        'PopulationSize', 100, 'Generations', 300);
    [x0, errfun] = ga(funGoal, nvars, A, b, [], [], lb, ub, [], [], optionsGA);
    
%      options = optimoptions('lsqnonlin',...
%         'Display','iter-detailed',...
%         ...%'PlotFcn',@optimplotfirstorderopt,...
%         'SpecifyObjectiveGradient',true,...
%         'CheckGradients',false,'FiniteDifferenceStepSize',1e-4,...
%         'Algorithm','trust-region-reflective',...%levenberg-marquardt
%         'OptimalityTolerance', 1e-14, 'FunctionTolerance', 1e-14,...
%         'MaxIterations', 20000, 'MaxFunctionEvaluations', 60000);
%     [k, cost] = lsqnonlin(funK, x0, lb, ub, options);
   
    disp('MultiSearch running...');
    opt = optimoptions('fminunc','Display', 'iter',...
       'FiniteDifferenceStepSize',1e-22, 'FiniteDifferenceType', 'central',...
       'Algorithm','trust-region', 'SpecifyObjectiveGradient',true,'HessianFcn','objective');
        % ','CheckGradients',true,
    %x = fmincon(funGoal,x0,A,b,[],[],lb,ub,[],opt)
    %Cproblem = createOptimProblem('fmincon', 'x0', x0, 'objective', funGoal, 'options', opt,...
    %    'Aineq', A, 'bineq', b, 'lb', lb, 'ub', ub);
    Cproblem = createOptimProblem('fminunc', 'x0', x0,'objective', funGoal,'options', opt);
    ms = MultiStart;
    [k,errormulti] = run(ms,Cproblem,200) %disp

%     [k, errormulti] = run(GlobalSearch('StartPointsToRun','bounds-ineqs',...
%         'FunctionTolerance', 1e-6, 'MaxWaitCycle', 20,...
%         'NumStageOnePoints', 200, 'NumTrialPoints', 1000, 'MaxTime', 60), Cproblem);
    end

    
    

%k - var;
function y = C_runge(k, T, C0, ssv, v, isExact, k2, k3)
    if(nargin <= 6)
        k2 = [];
        k3 = [];
    end
    if (~isExact)
        %Euler method
        size=length(T);
        y=zeros(1,size);
        for i=1:size
            if(i == 1)
                y(1)=C0;
            else
                y(i) = y(i-1) + F(0, y(i-1), k2, k3) * (T(i)-T(i-1));
            end
        end
    else
        %Runge Kutta
        options = odeset('RelTol',1e-4, 'AbsTol', 1e-4);
        [~,y] = ode23t(@(T, C)F(T, C, k2, k3), T, C0, options);
        y=y';
    end

function dCdt = F(t, c, k2, k3)
    if(isempty(k2) && isempty(k3))
        dCdt = -k(1)*ssv*c/((c+k(2)+c.^2/k(3)));
    else
        dCdt = -k(1)*ssv*c/((c+k2+c.^2/k3));
    end
end
end
    
function y = C_trap(k, T, C, ssv, v)
    size=length(T);
    y=zeros(1,size);
    for i=1:size
        if(i == 1)
            y(1)=C(1);
        else
            y(i)= y(i-1)+(F(k, i)+F(k, i-1))./2 .* (T(i)-T(i-1));
        end
    end
        function y = F(k, i)
            y = -k(1)*ssv*C(i)/((C(i)+k(2)+C(i).^2/k(3)));
        end
end

function data = augmentData(data, factor)
    C = data.C;
    T = data.T;
    X = data.X;
  
    [xData, yData] = prepareCurveData(T, X);
    ft = fittype('smoothingspline');
    opts = fitoptions('Method', 'SmoothingSpline', 'SmoothingParam', 0.9999);
    [fitresultX, gof] = fit(xData, yData, ft, opts);

    [xData, yData] = prepareCurveData(T, C);
    fitresultC = slmengine(xData,yData,'plot','on','decreasing','on', 'concavedown',...
        [0, max(C)*0.3],'rightminvalue',0); %
    %[fitresultC, gof] = fit(xData, yData, ft, opts);
    
    rate = 2;
%     for i = 1:factor
%         augSize=rate*length(T) -  1;
%         Taug = zeros(1, augSize);
%         Caug = zeros(1, augSize);
%         Xaug = zeros(1, augSize);
% 
%         Taug(1:rate:augSize) = T;
%         Caug(1:rate:augSize) = C;
%         Xaug(1:rate:augSize) = X;
% 
%         Taug(2:rate:augSize) = T(1:end-1)+diff(T)/2;
%         Caug(2:rate:augSize) = slmeval(Taug(2:rate:augSize),fitresultC);
%         Xaug(2:rate:augSize) = fitresultX(Taug(2:rate:augSize));
%         T = Taug;
%         C = Caug;
%         X = Xaug;
%     end
    augSize=rate.^factor*length(T) -  1;
    Taug = linspace(T(1),T(end), augSize);
    Caug = slmeval(Taug, fitresultC);
    Xaug = fitresultX(Taug);

    
    % Plot fit with data.
%     figure( 'Name', 'untitled fit 1' );
%     h = plot( fitresultC, Taug, Caug );
    hold on;
    plot( Taug, Caug, 'bo');
    %legend( h, 'C vs. T', 'untitled fit 1', 'Location', 'NorthEast' );
    % Label axes
    xlabel T
    ylabel C
    grid on
    data.C = Caug;
    data.X = Xaug;
    data.T = Taug;
    
end

%         indAug = round(length(T)/3);
%         Taug = [linspace(T(1), T(indAug), 10)'; T(indAug+1:end)];