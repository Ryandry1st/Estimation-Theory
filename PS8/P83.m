%% PS SPF #3
clear all;
close all;

plotloc1 = [0, 1600, 900, 600];
plotloc2 = [900, 1600, 900, 600];
plotloc3 = [1800, 1600, 900, 600];
plotloc4 = [0, 300, 900, 600];
plotloc5 = [900, 300, 900, 600];
plotloc6 = [1800, 300, 900, 600];

% xdot = s * cos(theta) + vx
% ydot = s * sin(theta) + vy
% sdot = vs
% thetadot = vtheta
% ldot = vl
% wdot = vw

% state = [x; y; s; theta; l; w]

plotting = 10;

dt = 0.1;
Q = diag([0.25, 0.25, 3, 40*pi/180, 0.1, 0.1])^2/dt;
tk = -dt;

% positive bearing is counter clockwise
% lidar is at origin
% z = [bmin, bmax, rmin]' + R
R = diag([2*pi/180, 2*pi/180, 0.1])^2;

load('problem3dataMod.mat');
options = odeset('reltol',1e-8);

nx = 6;
nv = 6;
nz = 3;

x(:, 1) = [90 4.25 13 pi 5 2]';
P(:, :, 1) = diag([2, 5, 1, pi/4, 4, 2])^2;

z = lidar.z;

beta = 2;
kappa = 0;
alpha = 1e-3;


lambda = alpha^2*(nx + nv + kappa) - (nx+nv);

%% Begin Iteration
for time=1:size(lidar, 1)
    Xa = [x(:, time); zeros(6, 1)];
    Pa = [P(:, :, time), zeros(6, 6); zeros(6, nx), Q];
    Sx = chol(Pa)';

    chi(:, 1) = Xa;
    for i=2:2*(nx+nv)+1
        if i <= nx+nv+1
            chi(:, i) = chi(:, 1) + sqrt(nx+nv+lambda)*Sx(:, i-1);
        else
            chi(:, i) = chi(:, 1) - sqrt(nx+nv+lambda)*Sx(:, i-nx-nv-1);
        end
    end

    chibar = zeros(12, 25);

    for k=1:25
        [TVEC, XMAT] = ode45(@dyn_car, [tk,tk+dt/2,tk+dt],chi(:, k),options);
        chibar(:, k) = XMAT(3, :);
    end

    w0m = lambda/(nx+nv+lambda);
    w0c = w0m + 1-alpha^2+beta;
    for i=2:2*(nx+nv)+1
        w(1, i) = 1/(2*(nx+nv+lambda));
        w(2, i) = w(1, i);
    end
    w(1, 1) = w0m;
    w(2, 1) = w0c;

    xbar = zeros(nx+nv, 1);
    Pbar = zeros(nx+nv, nx+nv);
    xbar = w0m*chibar(:, 1);
    for i=2:2*(nx+nv)+1
        xbar = xbar + (0.5/(nx+nv+lambda))*chibar(:, i);
    end
    Pbar = w0c*(chibar(:, 1)-xbar)*(chibar(:, 1)-xbar)';
    for i=2:2*(nx+nv)+1
        Pbar = Pbar + w(2, i)*(chibar(:, i) - xbar)*(chibar(:, i) - xbar)';
    end
    xbartrue = xbar(1:6);
    Pbartrue = Pbar(1:nx, 1:nx);
    %% measurement update

    xabar = [xbartrue; zeros(3, 1)];
    Pabar = [Pbartrue, zeros(nx, nz); zeros(nz, nx), R];

    Sxbar = chol(Pabar)';

    lambda = alpha^2*(nx+nz+k)-(nx+nz);

    g(:, 1) = xabar;
    for i=2:2*(nx+nz)+1
        if i <= nx+nz+1
            g(:, i) = g(:, 1) + sqrt(nx+nz+lambda)*Sxbar(:, i-1);
        else
            g(:, i) = g(:, 1) - sqrt(nx+nz+lambda)*Sxbar(:, i-nx-nz-1);
        end
    end

    tk = -dt;
    gbar = zeros(3, 2*(nx+nz)+1);

    for k=1:2*(nx+nz)+1
        gbar(:, k) = h_car(g(:, k));
    end

    w0m = lambda/(nx+nz+lambda);
    w0c = w0m + 1-alpha^2+beta;
    for i=2:2*(nx+nz)+1
        w(1, i) = 1/(2*(nx+nz+lambda));
        w(2, i) = w(1, i);
    end
    w(1, 1) = w0m;
    w(2, 1) = w0c;

    zbar = zeros(nz, 1);
    Pzz = zeros(nz, nz);

    zbar = w0m*gbar(:, 1);
    for i=2:2*(nx+nz)+1
        zbar = zbar + (0.5/(nx+nz+lambda))*gbar(:, i);
    end

    Pzz = w0c*(gbar(:, 1) - zbar)*(gbar(:, 1) - zbar)';
    for i=2:2*(nx+nz)+1
        Pzz = Pzz + w(2, i)*(gbar(:, i) - zbar)*(gbar(:, i) - zbar)';
    end

    Pxz = zeros(nx, nz);
    Pxz = w0c*(g(1:nx, 1) - xbartrue)*(gbar(:, 1) - zbar)';
    for i=2:2*(nx+nz)+1
        Pxz = Pxz + w(2, i)*(g(1:nx, i) - xbartrue)*(gbar(:, i) - zbar)';
    end

    %% Record new state and covariance
    zmeasure = [min(lidar(time).z(:, 2)), max(lidar(time).z(:, 2)), min(lidar(time).z(:, 1))]';
    x(:, time+1) = xbartrue + Pxz*inv(Pzz)*[zmeasure-zbar];
    P(:, :, time+1) = Pbartrue - Pxz*inv(Pzz)*Pxz';
end

load('problem3truth.mat')
for timestep=1:size(car, 1)
    xtruth(:, timestep) = car(timestep).x;
end


iterator = 1:220;
figure
hold on;
plot(iterator, xtruth(1, :), 'r', iterator, x(1, 2:end)', 'b');
plot(iterator, x(1, 2:end) + sqrt(squeeze(P(1, 1, 2:end)))', 'Color', [17 17 17]/256);
plot(iterator, x(1, 2:end) - sqrt(squeeze(P(1, 1, 2:end)))', 'Color', [17 17 17]/256);
title("Comparison of estimated and truth data for x position")
legend("Truth data", "Estimator");
set(gcf, 'position', plotloc1);


if plotting == 10
    figure
    hold on;
    plot(iterator, xtruth(2, :), 'r', iterator, x(2, 2:end)', 'b');
    legend("Truth data", "Estimator");
    plot(iterator, x(2, 2:end) + sqrt(squeeze(P(2, 2, 2:end)))', 'Color', [17 17 17]/256);
    plot(iterator, x(2, 2:end) - sqrt(squeeze(P(2, 2, 2:end)))', 'Color', [17 17 17]/256);
    legend("Truth data", "Estimator");
    title("Comparison of estimated and truth data for y position")
set(gcf, 'position', plotloc2);

    figure
    hold on;
    plot(iterator, xtruth(3, :), 'r', iterator, x(3, 2:end)', 'b');
    plot(iterator, x(3, 2:end) + sqrt(squeeze(P(3, 3, 2:end)))', 'Color', [17 17 17]/256);
    plot(iterator, x(3, 2:end) - sqrt(squeeze(P(3, 3, 2:end)))', 'Color', [17 17 17]/256);
    title("Comparison of estimated and truth data for speed")
    legend("Truth data", "Estimator");
set(gcf, 'position', plotloc3);

    figure
    hold on;
    x(4, 2:end) = mod(x(4, 2:end), 2*pi);
    plot(iterator, xtruth(4, :), 'r', iterator, x(4, 2:end)', 'b');
    plot(iterator, x(4, 2:end) + sqrt(squeeze(P(4, 4, 2:end)))', 'Color', [17 17 17]/256);
    plot(iterator, x(4, 2:end) - sqrt(squeeze(P(4, 4, 2:end)))', 'Color', [17 17 17]/256);
    title("Comparison of estimated and truth data for bearing")
    legend("Truth data", "Estimator");
set(gcf, 'position', plotloc4);

    figure
    hold on;
    plot(iterator, xtruth(5, :), 'r', iterator, x(5, 2:end)', 'b');
    plot(iterator, x(5, 2:end) + sqrt(squeeze(P(5, 5, 2:end)))', 'Color', [17 17 17]/256);
    plot(iterator, x(5, 2:end) - sqrt(squeeze(P(5, 5, 2:end)))', 'Color', [17 17 17]/256);
    title("Comparison of estimated and truth data for length")
    legend("Truth data", "Estimator");
    set(gcf, 'position', plotloc5);
    
    figure
    hold on;
    plot(iterator, xtruth(6, :), 'r', iterator, x(6, 2:end)', 'b');
    plot(iterator, x(6, 2:end) + sqrt(squeeze(P(6, 6, 2:end)))', 'Color', [17 17 17]/256);
    plot(iterator, x(6, 2:end) - sqrt(squeeze(P(6, 6, 2:end)))', 'Color', [17 17 17]/256);
    title("Comparison of estimated and truth data for width")
    legend("Truth data", "Estimator");
    set(gcf, 'position', plotloc6);
    
end

if plotting == 1
    h = figure;
    for i=1:220
        plotcar(x(:, i), '.-b', h)
        plotcar(xtruth(:, i), '.-r', h)
    end
end

% Large Alpha -- bad in general except width
% Small Alpha -- Normal
% Small Beta -- Speed is inverted, width is worse. bearing off by 180*
% Large Beta -- speed worse, things shift but generally same performance
% Large Kappa -- nothing changes, speed slightly worse
% Negative Kappa -- Nothing
