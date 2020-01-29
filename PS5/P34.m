% Problem 3 and 4 from PS 5
%% Solve A SLTI Kalman Filter problem
% Load data
clear all;
close all;

kf_example02a;
Qk = 10;
Rk = 0.025;

xhat = zeros(51, 2);
xhat(1, :) = xhat0;
P = zeros(51, 2, 2);
P(1, :, :) = P0;

ev = zeros(50);

%% Begin Filter
for k=1:50
    %% State and Cov. Prop
    xbar = Fk*xhat(k, :)';
    Pbar = Fk*squeeze(P(k, :, :))*Fk' + Gammak*Qk*Gammak';
    %% Measurement Update
    nu = zhist(k) - Hk*xbar;
    s = Hk*Pbar*Hk' + Rk;
    W = Pbar*Hk'*s^-1;
    xhat(k+1, :) = xbar + W*nu;
    P(k+1, :, :) = Pbar - W*s*W';
    ev(k) = nu*inv(s)*nu';
end

figure
plot(thist, xhat(2:end, 1), 'x', thist, sqrt(P(2:end, 1, 1)), 'o')
legend("Xhat(1)", "sqrt(P11(k))");
title("Prediction of X1 and prediction std")
xlabel("thist")

figure
plot(thist, xhat(2:end, 2), 'x', thist, sqrt(P(2:end, 2, 2)), 'o')
legend("Xhat(2)", "sqrt(P22(k))")
title("Prediction of X2 and prediction std")
xlabel("thist")

figure
plot(ev(:, 1), 'o')
legend("ev")
xlabel("k")
title("ev vs k")


P50 = squeeze(P(51, :, :))
Xhat50 = squeeze(xhat(51, :))

%% Question 4
% use kalman.m
% 
% sys = ss(Fk,[zeros(2, 1) Gammak], Hk, 0, -1);
% [km,L, Pss,Wss] = kalman(sys,Qk,Rk);
% Pss
% Wss
% error_trans = [eye(2) - Wss*Hk]*Fk;
% abs(eig(error_trans))