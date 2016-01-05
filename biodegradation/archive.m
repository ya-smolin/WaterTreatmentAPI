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