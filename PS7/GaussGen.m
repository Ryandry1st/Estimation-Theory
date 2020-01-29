function [x] = GaussGen(xbar,Pxx, nx, n)
% Generate nx, n gaussians with Pxx and xbar

A = chol(Pxx)';
%[ V, D] = eig(Pxx);
% A = real(V*sqrt(D));
x = A*randn(nx, n) + xbar;

% for i=1:n
%     z = randn(nx, 1);
%     % Use eigenvalue decomposition to solve for A for Pxx
%     
%     x(:, i) = A*z + xbar;
% end

