function [fxg] = p4_1(xg)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
x1 = xg(1);
x2 = xg(2);
top = x1 + x2 + x1*x2+5;
bot = x1^2+2*x2-x2^2-2;
fxg = [top; bot];
normal = norm(fxg, 2);
end

