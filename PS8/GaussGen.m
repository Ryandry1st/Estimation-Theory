function [x] = GaussGen(xbar,Pxx, nx, n)
% Generate nx, n gaussians with Pxx and xbar
A = chol(Pxx)';
x = A*randn(nx, n) + xbar;


