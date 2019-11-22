% Complete Problem 4-4
%% Initial setup
% model is z = h(x) + w
% Use bearing in the cost for comparison, but do not add information yet
clear all;
close all;

nk = 28;
nz = 2;
alpha = 0.2;
% 
% thist = [5;15;25;35;45;55;65;75;85;95;105;115];
% rhoahist = [352862.778204296;340399.972183495;328375.512355275;316754.551972577;305482.784782397;294500.455434987;283702.165996909;273085.433533899;262551.991714750;252050.375922058;241496.208747970;230858.792975043];
% rhobhist = [382864.360971427;370335.998477794;358259.255367332;346495.952598628;335082.681393646;323823.527400930;312866.241273828;301958.579383796;291075.842630371;280329.997342317;269465.468386008;258410.947432348];
load('radarmeasdata_missile_new.mat');
z = [rhoahist; rhobhist];

xg = P4_4_init(rhoahist, rhobhist, 3e4, thist);

[initial_cost, x, y, ra_hat, rb_hat] = Cost4_4(xg, thist, rhoahist, rhobhist, 1, 1, nk);
%% Iterate to find best xg

for i=1:50    
    [subh, ~ , ~, subH] = H_4_4(xg, thist, nk, nz); % H is 12x4x2
    
    % initial cost
    if i == 1
        Cost4_4(xg, thist, rhoahist, rhobhist, 1, 1, nk);
    end
    
    %% Update guess
    xg = xg + alpha * inv(subH'*subH)*subH'*[z-subh];
    
end

%% Results

[cost, x, y, ra_hat, rb_hat] = Cost4_4(xg, thist, rhoahist, rhobhist, 1, thetabhist, nk);

figure
plot(thist, rhoahist, thist, ra_hat);
title("Ra Comparison");
legend("Measured", "Predicted");

figure
plot(thist, rhobhist, thist, rb_hat);
title("Rb Comparison");
legend("Measured", "Predicted");

R = diag([10^2*ones(1,nk), 30^2*ones(1, nk)]);

format short
Pxx = inv(subH'*inv(R)*subH)
format long
cost
xg



%% Add bearing measurement
disp("Adding bearing information and recalculating")

nz = 4; % four usable measurements

xg = P4_4_init(rhoahist, rhobhist, 3e4, thist);
z = [rhoahist; rhobhist; thetaahist; thetabhist];

[initial_cost, ~, ~, ~, ~, ~, ~] = Cost4_4(xg, thist, rhoahist, rhobhist, thetaahist, thetabhist, nk);
%% Iterate to find best xg

for i=1:50    
    [subh, ~ , ~, subH] = H_4_4(xg, thist, nk, nz); % H is 12x4x2
    
    % initial cost
    if i == 1
        Cost4_4(xg, thist, rhoahist, rhobhist, thetaahist, thetabhist, nk);
    end
    
    %% Update guess
    xg = xg + alpha * inv(subH'*subH)*subH'*(z-subh);
    
end

[cost, x, y, ra_hat, rb_hat, ta_hat, tb_hat] = Cost4_4(xg, thist, rhoahist, rhobhist, thetaahist, thetabhist, nk);

figure
plot(thist, thetaahist, thist, ta_hat, '--');
title("Theta A Comparison");
legend("Measured", "Predicted");

figure
plot(thist, thetabhist, thist, tb_hat, '--');
title("Theta B Comparison");
legend("Measured", "Predicted");

R = diag([10^2*ones(1, nk), 30^2*ones(1, nk), 0.01^2*ones(1, nk), 0.03^2*ones(1, nk)]);

format short
Pxx = inv(subH'*inv(R)*subH)
format long
cost
xg


%% Bearing only measurement

% xg = xg*1.01;
xg = P4_4_init(rhoahist, rhobhist, 3e4, thist);
z = [thetaahist; thetabhist];

for i=1:50    
    [subh, ~ , ~, subH] = H_4_4ext(xg, thist, nk, nz); % H is 12x4x2
    
    % initial cost
    if i == 1
        Cost4_4(xg, thist, rhoahist, rhobhist, thetaahist, thetabhist, nk);
    end
    
    %% Update guess
    delta = inv(subH'*subH)*subH'*(z-subh);
    jold = Cost4_4(xg, thist, rhoahist, rhobhist, thetaahist, thetabhist, nk);
    jnew = Cost4_4(xg+alpha*delta, thist, rhoahist, rhobhist, thetaahist, thetabhist, nk);
    
    while jnew >= jold + 1e-3
        alpha = alpha/2.0;
        jnew = Cost4_4(xg+alpha*delta, thist, rhoahist, rhobhist, thetaahist, thetabhist, nk);
    end
    
    xg = xg + alpha * inv(subH'*subH)*subH'*(z-subh);
    
end

[cost, x, y, ra_hat, rb_hat, ta_hat, tb_hat] = Cost4_4(xg, thist, rhoahist, rhobhist, thetaahist, thetabhist, nk);

R = diag([0.01^2*ones(1, nk), 0.03^2*ones(1, nk)]);

format short
Pxx = inv(subH'*inv(R)*subH)
format long
cost
xg

rank(subH)