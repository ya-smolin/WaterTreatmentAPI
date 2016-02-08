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