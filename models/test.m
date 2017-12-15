clear
M = 88;
N = 100;

Tmax=1;
Hmax=1;
theta = Tmax / (N - 1);
h = Hmax / (M - 1);

y = zeros(M, N); %solution
y(:, 1) = 0;
y(1, :) = 1;
y(M, :) = 0;

beta = 1/h^2;
alpha = 1/theta;
nodesPerLayer = M-2;
upd=1;
diag = -2*(1+alpha/beta);
lowd = 1;
A = full(gallery('tridiag',nodesPerLayer, upd, diag, lowd));

for layer = 2:N
    b = zeros(nodesPerLayer, 1);
    for i = 2:M-1
        b(i-1) = 2*(1-alpha/beta)*y(i, layer-1) -y(i+1, layer-1) -y(i-1, layer-1);
        if i == 2
            b(i-1) = b(i-1) - y(1, layer);
        elseif i == M-1
            b(i-1) = b(i-1) - y(M, layer);
        end
    end
    yy = A\b;
    y(2:M-1, layer) = yy;
end

for x = 2:M
    %plot(0:h:Hmax,  y(:,layer), 'r-');
    plot(0:theta:Tmax,  y(x,:), 'r-');
    hold on;
    pause(0.5);
    cla;
end
    
function aa= superComplexSolution()
    M = 88;
    N = 100;

    Tmax=1;
    Hmax=1;
    theta = Tmax / (N - 1);
    h = Hmax / (M - 1);


    y = zeros(M, N); %solution
    y(:, 1) = 0;
    y(1, :) = 1;
    y(M, :) = 0;


    beta = 1/h^2;
    alpha = 1/(2*theta);

    %
    linInd=@(i,j)(j-2).*(M-2) + i-1; %linear index of the grid point
    linDim = linInd(M-1, N); %last node lin index
    isNode = @(i,j)linInd(i,j) > 0 && linInd(i,j) <= linDim;
    A = zeros(linDim, linDim);%linear system for noders
    b = zeros(linDim, 1);
    for j = 2:N
        for i = 2:M-1
            eqNum = linInd(i,j);
            b(eqNum) = 0;
            if(j == N)
               processNode(i, j, -2*beta-2*alpha);
               processNode(i, j-1, 2*alpha);
            else
               processNode(i, j, -2*beta);
               processNode(i, j-1, alpha);
               processNode(i, j+1, -alpha);
            end    
            processNode(i-1, j, beta);
            processNode(i+1, j, beta);
        end
    end
    yy = A \ b;
    y(2:M-1, 2:N) = reshape(yy, M-2, N-1);

    for i = 1:M
        plot(0:theta:Tmax,  y(i,:), 'r-');
        pause(0.05);
        if i ~= M
            cla;
        end
    end
    for j = 1:N
        plot(0:h:Hmax,  y(:,j), 'r-');
        pause(0.5);
        if j ~= N
            cla;
        end
    end

    function processNode(i, j, val)
    %     disp([num2str(i) ' ' num2str(j)]);
        if isNode(i,j)
            A(eqNum, linInd(i,j)) = val;
        else
            b(eqNum) = b(eqNum) - val * y(i,j);
        end
    end
end

function t=page315318
    n=10;
    h=1/n;
    y0 = 1;
    alpha=1/h^2+1/(2*h);
    beta = -10-2/h^2;
    gamma = 1/h^2-1/h;
    A = full(gallery('tridiag',n, alpha, beta, gamma));
    A(1,1:2) = [beta, gamma];
    A(n, n-1:n) = [-1/h 1/h];
    b = zeros(n, 1);
    b(1) = -y0*alpha;
    y = A\b
    plot(linspace(0,1,n+1), [y0; y], 'r-')
end


function t = test1
Tmax = 100;
N = 100;
y = zeros(1,N+1);

%initial condtions
y0=1;
y(1)=y0;
h = Tmax/N;

%Euler explicit
for i = 2:N+1
    y(i) = y(i-1)-h*(y(i-1)^2);
end

hold on;
fplot(@(t)(1./(t+1)), [0, Tmax],'DisplayName',['exact solution, h=' num2str(h)]);
plot(0:h:Tmax, y, 'bo','DisplayName','Euler explicit');
y = fsolve(@(y)modelEuler(y, y0, N, h), ones(1,N));
plot(0:h:Tmax, [y0 y], 'r+','DisplayName','Euler implicit');
y = fsolve(@(y)modelTrap(y, y0, N, h), ones(1,N));
plot(0:h:Tmax, [y0 y], 'g*','DisplayName','Trap');
legend();
end

%Euler implicit
function F = modelEuler(y, y0, N, h)
    F = zeros(N, 1);
    for i = 1:N
        if(i==1)
            F(1) = h*y(1).^2 + y(1) - y0; 
        else
            F(i) = h*y(i).^2 + y(i) - y(i-1);
        end
    end
end

%Euler implicit
function F = modelTrap(y, y0, N, h)
    F = zeros(N, 1);
    for i = 1:N
        if(i==1)
            F(1) = (y(1)-y0)./h+((y0+y(1))./2).^2;
        else
            F(i) = (y(i)-y(i-1))./h+((y(i-1)+y(i))./2).^2;
        end
    end
end