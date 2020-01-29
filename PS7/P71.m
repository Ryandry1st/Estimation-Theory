%% Solving 7.1
clear all;
format long;

uk = zeros(3, 1);
vk = zeros(3, 1);
% norm(rs) = 6478+780;
dist = 6378e3+780e3;


mu = 3.986005e14;
% rx = dist;
% ry = 3000;
% rz = sqrt(dist^2 - rx^2 - ry^2);
% rs = [rx, ry, rz]';
% vs = [-ry, rx, 0]';

rx = dist;
ry = 0;
rz = 0;
rs = [rx, ry, rz]';
vs = [0, sqrt(mu/dist), 0]';

x = [rs; vs];

tk = 0;

% total time = distance/speed
% tt = round(2*pi*sqrt(dist^3/mu));
tt = 2048;
xhat = x;
T = 1;
T= 2*pi*dist/(norm(vs)*tt);

x = zeros(6, round(tt));
x(:, 1) = xhat;
for i=1:tt
    [x(:, i+1), fk, gk] = propagateOrbit(tk,T, x(:, i), uk, vk, mu);
    diff(i) = norm(x(:, i+1)-xhat);
    tk = tk + T;
end

 change = xhat - x(:, end)
 norm(change)
 
 % plot(diff)
 [minval, minarg] = min(diff);
 
 plot3(x(1, :), x(2, :), x(3, :)); hold on
 [xdumb, ydumb, zdumb] = sphere; surf(xdumb.*6378e3, ydumb.*6378e3, zdumb.*6378e3);

%% Test F
xbar = xhat*.999;
T = 1;
[xhat2, fk2, gk2] = propagateOrbit(tk,T, xhat, uk, vk, mu);
[xbar2, fk3, gk3] = propagateOrbit(tk,T, xbar, uk, vk, mu);
lhs = [xhat2 - xbar2];
rhs = fk2*[xhat - xbar];
lhs - rhs
