function [xg_new] = gradient_descent(xg)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
x1 = xg(1);
x2 = xg(2);
V = [1 + x2, 1+x1; 2*x1, 2-2*x2];
xg_new = xg - inv(V)*p4_1(xg);

end

