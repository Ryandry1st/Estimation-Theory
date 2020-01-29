function [xk] = Propagate(xkm1, ukm1, vkm1)
x = xkm1(1);
y = xkm1(2);
theta = xkm1(3);

ds = ukm1(1) + vkm1(1);
dtheta = ukm1(2) + vkm1(2);

% first rotate, then move
xk = zeros(3, 1);
xk(3) = theta + dtheta;
xk(1) = x + ds*cos(xk(3));
xk(2) = y + ds*sin(xk(3));
end

