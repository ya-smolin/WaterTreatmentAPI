%clear;

v=2.1;
ssv = 825;
k = [4.65 17.1 31.4];
T = linspace(0, 350/(60*24), 20)
res = @(t, y)-k(1)*ssv*y/(v*(y+k(2)+y.^2/k(3))); %Haldane

options = odeset('RelTol',1e-4,'AbsTol',[1e-4]);
[tt,Y] = ode45(res,T,[132],options);
plot(tt, Y, 'r-');