% Problem Set 6 number 3
% Smoothing filter
% Calculate the smoothed estimates for the problem in kf example03a.m. 
% Compare x?(10) with x?(10) and compare P (10) with P ?(10). 
% Is P ?(10) ? P (10)? Do the smoothed state time history estimate plots
% look ?smoother? than the filtered state time history estimate plots?

clear all;
close all;

kf_example03a;

plotting1 = 0;
plotting2 = 1;

xhat = zeros(51, 3);
%% Begin SRIF
returnx = 1;

% Add g and u for versatility
G = 0;
u = zeros(3, 1);

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
wx(1, :) = zx(1, :)' - (squeeze(Rxx(1, :, :))*xhat0);
xhat2(1, :) = inv(squeeze(Rxx(1, :, :)))*zx(1, :)';

% static Rk
Ra = chol(Rk);
invRat = inv(Ra)';

for k=1:50
    % propagation step
    [Qa, RA] = qr([Rvv, zeros(nv, nx); -1*squeeze(Rxx(k, :, :))*inv(Fk)*Gammak, squeeze(Rxx(k, :, :))*inv(Fk)]);
    dumbmatrixtoseperate = Qa'*[zeros(nv, 1); zx(k, :)'+ (squeeze(Rxx(k, :, :))*G*u)];
    zvbar(k, :) = dumbmatrixtoseperate(1:nv);
    zxbar(k, :) = dumbmatrixtoseperate(nv+1:end);
    Rxxbar(k, :, :) = squeeze(RA(nv+1:end, nv+1:end));
    Rvxbar(k, :, :) = squeeze(RA(1:nv, nv+1:end));
    Rvvbar(k, :, :) = squeeze(RA(1:nv, 1:nv));
    % measurement update
    za = invRat*zhist(k);
    Ha = invRat*Hk;
    
    [Qb, Rb] = qr([squeeze(Rxxbar(k, :, :)); Ha]);
    dumbmatrix2 = Qb'*[zxbar(k, :)'; za'];
    zx(k+1, :) = dumbmatrix2(1:nx);
    Rxx(k+1, :, :) = Rb(1:nx, :);
    if returnx == 1
        xhat2(k+1, :) = inv(squeeze(Rxx(k+1, :, :)))*zx(k+1, :)';
        Pxx2(k+1, :, :) = inv(squeeze(Rxx(k+1, :, :)))*inv(squeeze(Rxx(k+1, :, :)))';
        v(k, :) = xhat2(k+1, :)' - Fk*xhat2(k, :)' - G*u;
        % wx(k+1, :) = zx(k+1, :)'-(squeeze(Rxx(k+1, :, :))*xhat2(k+1, :)');
    end
    % wvbar(k, :) = zvbar(k, :) - squeeze(Rvvbar(k, :, :))*v(k)
end
xhat50_srif = xhat2(end, :)
P50_srif = squeeze(Pxx2(end, :, :)) 

%% plot the results
if plotting1==1
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

%% Smoothing Filter
zxstar = zeros(size(zx));
Rxxstar = zeros(size(Rxx));
Wxstar = zeros(size(zx));
xstar = zeros(size(xhat2));
Pstar = zeros(size(Pxx2));
zvstar = zeros(50, nv);

zxstar(end, :) = zx(end, :);
Rxxstar(end, :, :) = Rxx(end, :, :);
% Wxstar(end, :) = (zx(end, :)' - squeeze(Rxx(end, :, :))*xhat2(end, :, :)')';

xstar(end, :) = inv(squeeze(Rxxstar(end, :, :)))*zxstar(end, :)';
Pstar(end, :, :) = inv(squeeze(Rxxstar(end, :, :)))*inv(squeeze(Rxxstar(end, :, :)))';

% x(k+1) = F(k)*x(k) + G(k)*u(k) + L(k)v(k)
top = 50;

for k=top:-1:1
    if k<=top
        block1 = [squeeze(Rvvbar(k, :, :)) + squeeze(Rvxbar(k, :, :))*Gammak, squeeze(Rvxbar(k, :, :))*Fk; squeeze(Rxxstar(k+1, :, :))*Gammak, squeeze(Rxxstar(k+1, :, :))*Fk];
        [Qc, Rc] = qr(block1);
        rhs = Qc'* block1;
        lhs = Qc' * [zvbar(k, :)' - (squeeze(Rvxbar(k, :, :))*G*u); zxstar(k+1, :)' - (squeeze(Rxxstar(k+1, :, :))*G*u)];
        
    else
        block1 = [squeeze(Rvvbar(k, :, :)), squeeze(Rvxbar(k, :, :)); zeros(nx, nv), squeeze(Rxxstar(k+1, :, :))];
        [Qc, Rc] = qr(block1);
        rhs = Qc'* block1;
        lhs = Qc' * [zvbar(k, :)'; zxstar(k+1, :)'];
    end
    
    zvstar(k, :) = lhs(1:nv);
    zxstar(k, :) = lhs(nv+1:end);
    Rvvstar = rhs(1:nv, 1:nv);
    Rvxstar = rhs(1:nv, nv+1:end);
    Rxxstar(k, :, :) = rhs(nv+1:end, nv+1:end);
    xstar(k, :) = inv(squeeze(Rxxstar(k, :, :)))*zxstar(k, :)';
    vstar(k, :) = inv(Rvvstar)*(zvstar(k, :)' - Rvxstar*xstar(k, :)');
    Pstar(k, :, :) = inv(squeeze(Rxxstar(k, :, :)))*inv(squeeze(Rxxstar(k, :, :)))';
    
end

if plotting2 == 1
    figure
    plot(1:51, xhat2(:, 1), 1:51, xstar(:, 1))
    legend("Xhat", "Xstar");
    title("X(1)")

    figure
    plot(1:51, xhat2(:, 2), 1:51, xstar(:, 2))
    legend("Xhat", "Xstar");
    title("X(2)");

    figure
    plot(1:51, xhat2(:, 3), 1:51, xstar(:, 3))
    legend("Xhat", "Xstar");
    title("X(3)");


    figure
    plot(1:51, xhat2(:, 1) - xstar(:, 1), 1:51, xhat2(:, 2) - xstar(:, 2), 1:51, xhat2(:, 3) - xstar(:, 3));
    legend("X(1)", "X(2)", "X(3)");
    title("Difference of smoother vs estimator")

    x10 = xhat2(11, :) - xstar(11, :)
    P10 = squeeze(Pxx2(11, :, :) - Pstar(11, :, :))
    % if cholesky works it is positive definite
    chol(P10);
end
