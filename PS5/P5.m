% Problem set 5 # 5
% load the data
clear all
kf_example02b;
plotting = 0;

epsilon = zeros(3, 50);

xhat = zeros(3, 51, 2);
P = zeros(3, 51, 2, 2);
for j=1:3
    xhat(j, 1, :) = xhat0;
    P(j, 1, :, :) = P0;
end

for index=1:3

    % Begin Filter
    for k=1:50
        % State and Cov. Prop
        xbar = Fk*reshape(xhat(index, k, :), 1, 2)';
        Pbar = Fk*squeeze(P(index, k, :, :))*Fk' + Gammak*Qk(index)*Gammak';
        % Measurement Update
        nu = zhist(k) - Hk*xbar;
        s = Hk*Pbar*Hk' + Rk;
        W = Pbar*Hk'*s^-1;
        xhat(index, k+1, :) = xbar + W*nu;
        P(index, k+1, :, :) = Pbar - W*s*W';
        
        epsilon(index, k) = nu'*inv(s)*nu;
    end

    if plotting == 1
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
    end
   
end

%% calculate bounds for chi-square distribution for correct Q
alpha = 0.01;
r1 = chi2inv(alpha/2, 50)/50;
r2 = chi2inv(1-alpha/2, 50)/50;
ideal = 0;

means = mean(epsilon');
for j=1:3
    if means(j) >= r1
        if means(j) <= r2
            disp("The correct value for Q is ")
            disp(Qk(j))
            disp("which resulted in mean epsilon of ")
            disp(means(j))
            ideal = j;
        end
    end
end
switch ideal
    case 1
        notideal1 = 2;
        notideal2 = 3;
    case 2
        notideal1 = 1;
        notideal2 = 3;
    case 3
        notideal1 = 1;
        notideal2 = 2;
end


%% State estimate differences
differences = zeros(2, 40, 2);
differences(1, :, :) = xhat(ideal, 12:end, :) - xhat(notideal1, 12:end, :);
differences(2, :, :) = xhat(ideal, 12:end, :) - xhat(notideal2, 12:end, :);

% figure
% plot(thist, differences(1, :, :))

% figure

rms12 = sqrt(trace(squeeze(differences(1, :, :))'*squeeze(differences(1, :, :)))/40);
sigma12 = sqrt(trace(squeeze(P(ideal, 51, :, :))));
result(1) = rms12/sigma12;

rms32 = sqrt(trace(squeeze(differences(2, :, :))'*squeeze(differences(2, :, :)))/40);
sigma32 = sqrt(trace(squeeze(P(ideal, 51, :, :))));
result(2) = rms32/sigma32;
% Result(1) is less than 1 std away so it is fine
% Result(2) is 2.24 stds away so it is significant.