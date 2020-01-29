% Problem set 6 #1
% Develop a truth model MC simulation
clear all;
format long

MC = [50, 1000];
kmax = 50;
important1 = 10;
important2 = 35;

% load the basic params from Q3
kf_example02b
clear zhist
rng(0);

xhat = zeros(51, 2);
xhat(1, :) = xhat0;
P = zeros(51, 2, 2);
P(1, :, :) = P0;
xtilde = zeros(50, 2);
epsilon = zeros(3, 50, 1);
xplots = zeros(2, 1000, 50, 2);

%% Run Monte Carlos
for sets=1:2
    for sim=1:MC(sets)
        [xhist, zhist] = mcltisim(Fk, Gammak, Hk, Qk1, Rk, xhat0, P0, kmax);
        xhat(1, :) = xhat0;
        P(1, :, :) = P0;
        %% Begin Filter
        for k=1:kmax
            %% State and Cov. Prop
            xbar = Fk*xhat(k, :)';
            Pbar = Fk*squeeze(P(k, :, :))*Fk' + Gammak*Qk3*Gammak';
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
        xplots(sets, sim, :, :) = xtilde(:, :);
        xplots(sets, sim, :, :) = xtilde(:, :);
    end
%     disp(['MC=' num2str(MC(sets))]);
%     disp(['xtilde(10)=' num2str(xtilde(10))]);
%     disp(['xtilde(10)*xtildeT(10)=' num2str(xtilde(10)*xtilde(10)')]);
%     disp(['xtilde(35)=' num2str(xtilde(35))]);
%     disp(['xtilde*xtildeT(35)=' num2str(xtilde(35)*xtilde(35)')]);
end

figure
plot(1:kmax, xtilde)
title("Xtilde")

% legend(['MC=' num2str(MC(1))], ['MC=' num2str(MC(2))], ['MC=' num2str(MC(3))])
% plot(xhist)

special1 = mean(squeeze(xplots(2, :, important1, :)))
special2 = mean(squeeze(xplots(2, :, important2, :)))
cov1 = cov(squeeze(xplots(1, :, important1, :)))
disp("goes to ->")
cov1 = cov(squeeze(xplots(2, :, important1, :)))
disp("compare with")
squeeze(P(important1+1, :, :))
cov2 = cov(squeeze(xplots(1, :, important2, :)))
disp("goes to ->")
cov2 = cov(squeeze(xplots(2, :, important2, :)))
disp("compare with")
squeeze(P(important2+1, :, :))

figure
plot(squeeze(xplots(2, :, important1, :)))

figure
plot(squeeze(xplots(2, :, important2, :)))
