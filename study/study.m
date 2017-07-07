a = 45;
R = rotx(a);
Rp = R(2:3, 2:3);
R = rotx(-a);
Rm = R(2:3, 2:3);

N = 100;
xy = zeros(N, 2);
xy(1, :) = [1 0];
sum = xy;
for i = 2:N
    r = rand(1)
    xp = xy(i-1, :)';
    if r > 0.5
        Rx = Rp;
    else
        Rx = Rm;
    end
    xy(i, :) =   Rm * xp;
    sum(i, :) = xp +  xy(i, :)';
    arrow(sum(i-1, :), sum(i, :));
end
xx = hypot(sum(:,1), sum(:, 2));
 arrow([0 0], sum(1, :));

% t = linspace(0, N, N);
% plot(t, xx);