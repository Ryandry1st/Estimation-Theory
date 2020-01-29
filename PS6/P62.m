% Problem 2 from Problem set 6
% First KF then SRIF
clear all
close all
plotting = 1;

for example=1:2
    if example==1
        kf_example03a;
    else
        kf_example03b;
    end
    %% Solve with regular KF
    xhat = zeros(51, 3);
    xhat(1, :) = xhat0;
    P = zeros(51, 3, 3);
    P(1, :, :) = P0;

    %% Begin regular KF
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
    end
    xhat50 = xhat(end, :)
    P50 = squeeze(P(end, :, :))

    %% Begin SRIF
    returnx = 1;

    nv = size(Gammak, 2);
    nx = size(xhat, 2);

    xhat2 = zeros(51, nx);
    Rxx = zeros(51, nx, nx);
    Rvv = zeros(nv, nv);
    zx = zeros(51, nx);
    Pxx2 = zeros(51, 3, 3);


    Rxx(1, :, :) = [inv(chol(P0))]';
    Rvv = [inv(chol(Qk))]';
    zv = 0;
    zx(1, :) = squeeze(Rxx(1, :, :))*xhat0;

    % static Rk
    Ra = chol(Rk);
    invRat = inv(Ra)';

    for k=1:50
        % propagation step
        [Qa, RA] = qr([Rvv, zeros(nv, nx); -1*squeeze(Rxx(k, :, :))*inv(Fk)*Gammak, squeeze(Rxx(k, :, :))*inv(Fk)]);
        dumbmatrixtoseperate = Qa'*[zeros(nv, 1); zx(k, :)'];
        zxbar = dumbmatrixtoseperate(nv+1:end);
        Rxxbar = squeeze(RA(nv+1:end, nv+1:end, :, :));
        % measurement update
        za = invRat*zhist(k);
        Ha = invRat*Hk;
        [Qb, Rb] = qr([Rxxbar; Ha]);
        dumbmatrix2 = Qb'*[zxbar; za'];
        zx(k+1, :) = dumbmatrix2(1:nx);
        Rxx(k+1, :, :) = Rb(1:nx, :);
        if returnx == 1
            xhat2(k+1, :) = inv(squeeze(Rxx(k+1, :, :)))*zx(k+1, :)';
            Pxx2(k+1, :, :) = inv(squeeze(Rxx(k+1, :, :)))*inv(squeeze(Rxx(k+1, :, :)))';
        end
    end
    xhat50_srif = xhat2(end, :)
    P50_srif = squeeze(Pxx2(end, :, :)) 

    %% plot the results
    if plotting==1
        figure
        plot(thist, xhat(2:end, 1), 'x', thist, xhat2(2:end, 1), 'o')
        legend("KF", "SRIF");
        title("Prediction of X1 and prediction std")
        xlabel("thist")

        figure
        plot(thist, xhat(2:end, 2), 'x', thist, xhat2(2:end, 2), 'o')
        legend("KF", "SRIF")
        title("Prediction of X2 and prediction std")
        xlabel("thist")

        figure
        plot(thist, xhat(2:end, 3), 'x', thist, xhat2(2:end, 3), 'o')
        legend("KF", "SRIF");
        title("Prediction of X3 and prediction std")
        xlabel("thist")
    end
end
