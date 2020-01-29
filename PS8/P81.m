%% Problem Set 8 #1
clear all;
close all;

% three operating modes [c, hz]
op1 = [0.1, 1];
op2 = [0.5, 2];
op3 = [1, 3];
op = [op1; op2; op3];
sel = 1;

% xdot = [thetazdot, wzdot]' = [0, wz; 0, -c/hz*wz] + [0; 1/hz]*u
% z(tk) = [1 0]x(tk) + w(tk)
% x(tk) = [thetaz(tk), wz(tk)]' and w is zero mean, white w(tk)~N(0, R)

% Develop mulitple model filter to do estimation
% assume control input is zero order hold
% discretize dynamics to get x(k+1) = F(k)*x(k) + G(k)*u(k) + v(k)
A = [0, 1; 0, -op(sel, 1)/op(sel, 2)];

% A is constant so F is just expm(A), and G is just F*B
F = expm(A);
G = F*[0; 1/op(sel, 2)];
R = 0.1
Q = 0.001*diag([0.1, 1]);
x1 = [0; 0.1];


