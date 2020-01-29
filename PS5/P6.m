% Problem set 5 #6
% Develop a truth model MC simulation
clear all;
close all;

format short

% 5000 is used for the expectation case
maxMC = 5000;

MC = [50, 1000, maxMC];
kmax = 50;
time10 = 10;
time35 = 35;
% load the basic params from Q3
kf_example02a
clear zhist
rng(0);

xhat = zeros(51, 2);
xhat(1, :) = xhat0;
P = zeros(51, 2, 2);
P(1, :, :) = P0;
xtilde = zeros(50, 2);
epsilon = zeros(3, 50, 1);
xplots = zeros(1000, 50, 2);

%% Run Monte Carlos
for sets=1:3
    for sim=1:MC(sets)
        [xhist, zhist] = mcltisim(Fk, Gammak, Hk, Qk, Rk, xhat0, P0, kmax);
        xhat(1, :) = xhat0;
        P(1, :, :) = P0;
        %% Begin Filter
        for k=1:kmax
            %% State and Cov. Prop
            xbar = Fk*xhat(k, :)';
            Pbar = Fk*squeeze(P(k, :, :))*Fk' + Gammak*Qk*Gammak';
            %% Measurement Update
            nu = zhist(k) - Hk*xbar;
            s = Hk*Pbar*Hk' + Rk;
            W = Pbar*Hk'*s^-1;
            xhat(k+1, :) = xbar + W*nu;
            P(k+1, :, :) = Pbar - W*s*W';
            xtilde(k, :) = xhist(k+1, :) - xhat(k+1, :);
            % need to scale the epsilons by the number of MC run
            epsilon(sets, k) = (xtilde(k, :)*squeeze(P(k, :, :)) * xtilde(k, :)')./MC(sets);
        end 
        xplots(sim, :, :) = xtilde(:, :);
        xplots(sim, :, :) = xtilde(:, :);
    end
    disp(['MC=' num2str(MC(sets))]);
    
    disp(['E[xtilde(10)]=']);
    Extilde10 = mean(squeeze(xplots(:, time10, :)))
    
    disp(["E[xtilde(10)*xtilde'(10)]="])
    cov10 = cov(squeeze(xplots(:, time10, :)))
    disp("compare with P(10)")
    P10 = squeeze(P(time10+1, :, :))
    
    disp(['Extilde(35)]=']);
    Extilde35 = mean(squeeze(xplots(:, time35, :)))
    disp(["E[xtilde(35)*xtilde'(35)]="]);
    cov35 = cov(squeeze(xplots(:, time35, :)))
    disp("compare with P(35)")
    P35 = squeeze(P(time35+1, :, :))
end

figure
plot(1:kmax, xtilde)
title("Xtilde")
% 
% Extilde10 = mean(squeeze(xplots(:, time10, :)))
% Extilde35 = mean(squeeze(xplots(:, time35, :)))
% 
% cov10 = cov(squeeze(xplots(:, time10, :)))
% disp("compare with P(10)")
% P10 = squeeze(P(time10+1, :, :))
% 
% cov35 = cov(squeeze(xplots(:, time35, :)))
% disp("compare with P(35)")
% P35 = squeeze(P(time35+1, :, :))

sims = 1:maxMC;
figure

ave10 = zeros(maxMC, 2);
ave35 = zeros(maxMC, 2);
for i=1:maxMC
    ave10(i, :) = mean(squeeze(xplots(1:i, time10, :)));
    ave35(i, :) = mean(squeeze(xplots(1:i, time35, :)));
end

plot(sims, ave10)
title("Xtilde(10) average versus monte carlo runs");

figure
plot(sims, ave35)
title("Xtilde(35) average versus monte carlo runs");
