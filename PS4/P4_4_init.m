function [initial_cond] = P4_4_init(rhoa,rhob, dl, thist)
% assumes starting point for radar a is rad_a_x
rad_a_x = 4.1e5;
rad_b_x = rad_a_x + dl;

% General Method
theta = acos((-rhob.^2+rhoa.^2+dl^2)./(2.*rhoa*dl));
y = real(rhoa.*sin(theta));
x = real(rad_a_x + rhoa.*cos(theta));

%% select x(2), x(11), y(2), y(11) for parabola creation
x1 = x(2);
x2 = x(end/2+1);
h1 = y(2);
h2 = y(end/2+1);

a = -(h1-h2)/(x1 - x2)^2;
% the plot should be similar/fit on those points
% with parabola -a*(x-x2)^2+h2

figure
plot(x, y, x, -a.*(x-x2).^2 + h2);
title("X vs Y initial estimate")
legend("X,Y from two points", "X, Y parabola function");

%% initial condition when t = 0
% x is just initial position + initial velocity * time
delta_x = x2 - x1;
delta_t = thist(11) - thist(2);

Vx = delta_x/delta_t;
X_0 = x1 - Vx*thist(2);
Y_0 = -a*(X_0 - x2)^2+h2;
Vy = (h2 + 1/2*9.81*thist(11).^2-Y_0)/thist(11);

if Y_0 < 0
    Y_0 = 0;
end

initial_cond = [X_0, Y_0, Vx, Vy]';

end

