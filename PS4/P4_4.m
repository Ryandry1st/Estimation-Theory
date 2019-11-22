% Complete Problem 4-4
%% Initial setup
% model is z = h(x) + w
clear all;

nz = 12;
nx = 4;
alpha = 0.2;

thist = [5;15;25;35;45;55;65;75;85;95;105;115];
rhoahist = [352862.778204296;340399.972183495;328375.512355275;316754.551972577;305482.784782397;294500.455434987;283702.165996909;273085.433533899;262551.991714750;252050.375922058;241496.208747970;230858.792975043];
rhobhist = [382864.360971427;370335.998477794;358259.255367332;346495.952598628;335082.681393646;323823.527400930;312866.241273828;301958.579383796;291075.842630371;280329.997342317;269465.468386008;258410.947432348];
% load('radarmeasdata_missile_new.mat');
z = [rhoahist; rhobhist];

xg = P4_4_init(rhoahist, rhobhist, 3e4, thist);

[initial_cost, x, y, ra_hat, rb_hat] = Cost4_4(xg, thist, rhoahist, rhobhist);
%% Iterate to find best xg

for i=1:50
    HTH = zeros(4, 4);
    HTzh = zeros(4, 1);
    
    [subh, h , H, subH] = H_4_4(xg, thist); % H is 12x4x2
    
    % initial cost
    if i == 1
        Cost4_4(xg, thist, rhoahist, rhobhist);
    end
    zh = [z-h];
    
    %% define tensor multiplication
    % output of H'H = 4x4 and output of inv(H'H)*H' is 4x2x12
    % output of inv(H'H)*H'*[z-h] = 4x1
    for l=1:12
        HTH = HTH + squeeze(H(l, :, :))*squeeze(H(l, :, :))';
    end
    for l=1:12
        HTzh = HTzh + squeeze(H(l, :, :))* zh(l, :)';
    end
    
    %% Update guess
    xg = xg + alpha * inv(subH'*subH)*subH'*[z-subh];
    
end

[h, H] = H_4_4(xg, thist);
zh = [z-h];

%% Results

[cost, x, y, ra_hat, rb_hat] = Cost4_4(xg, thist, rhoahist, rhobhist);

figure
plot(thist, rhoahist, thist, ra_hat);
title("Ra Comparison");
legend("Measured", "Predicted");

figure
plot(thist, rhobhist, thist, rb_hat);
title("Rb Comparison");
legend("Measured", "Predicted");
R = diag([10^2*ones(1,12), 30^2*ones(1, 12)]);

Pxx = inv(subH'*inv(R)*subH);