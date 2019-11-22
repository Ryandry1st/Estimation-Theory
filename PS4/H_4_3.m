function [H] = H_4_3(xg, t)
x1 = xg(1);
x2 = xg(2);
x3 = xg(3);

H = zeros(3, 11);
for j=1:11
    H(:, j) = [cos(x2*t(j)+x3), -x1*t(j)*sin(x2*t(j)+x3), -x1*sin(x2*t(j)+x3)];
end

