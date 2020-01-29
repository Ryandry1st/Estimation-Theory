%% PS SPF 4
% Bootstrap Particle Filter

clear all;
close all;

plotloc1 = [0, 1600, 900, 600];
plotloc2 = [900, 1600, 900, 600];
plotloc3 = [1800, 1600, 900, 600];
plotloc4 = [0, 300, 900, 600];
plotloc5 = [900, 300, 900, 600];
plotloc6 = [1800, 300, 900, 600];

plotting = 1;

dt = 0.1;
% noise for how far moved and how far turned
Q = diag([0.1, 5*pi/180])^2;

% turns, moves, listens to sonar, doesnt know which beacons are which
R = diag([1, 1, 1])^2;

% state is [x, y, theta]'
% theta is measured counter-clockwise from the positive x direction
nx = 3;
nu = 2;
nz = 3;
num_beacons = 5;

particles = 1000;
timesteps = 101;
part_set = zeros(3, particles);
old_parts = zeros(timesteps, 3, particles);

mhist = zeros(timesteps, particles);
% prop_part_set = zeros(size(part_set));
resample_index = 500;

xhat = zeros(nx, timesteps);
P = zeros(nx, nx, timesteps);

load('problem4data.mat');

xset = minx + (maxx-minx).*rand(particles, 1);
yset = miny+ (maxy-miny).*rand(particles, 1);

thetaset = 2*pi.*rand(particles, 1);

% equal weighting
weights = ones(particles, 1)./particles;
logwiprimek = zeros(particles, 1);

for k=1:particles
    part_set(:, k) = [xset(k); yset(k); thetaset(k)];
end

for k=1:timesteps
    uk = encoder(k).u;
    vkminus1 = GaussGen(0, Q, nu, particles);
    % propagate the particles forward

    old_parts(k, :, :) = part_set(:, :);
    
    for j=1:particles
        part_set(:, j) = Propagate(part_set(:, j), uk, squeeze(vkminus1(:, j)));
    end

    
    zk = sonar(k).z;
    % compute ranges to the three closest beacons for each particle
    part_dist = zeros(num_beacons, particles);
    % initialize so that it can be overwritten easily
    nzbeacons = zeros(nz, particles);
    
    for j=1:particles
        for i=1:num_beacons
            % use euclidean distance i.e. norm
            part_dist(i, j) = norm(part_set(1:2, j) - beacons(i, :)');
        end
        nzbeacons(:, j) = mink(part_dist(:, j), nz);
        logwiprimek(j) = -1/2*(zk - nzbeacons(:, j))'*inv(R)*(zk - nzbeacons(:, j)) + log(weights(j));
    end
    
    % update weights via log likelihood and max
    maxlogwik = max(logwiprimek);
    weights(:) = exp(logwiprimek(:) - maxlogwik);
    weights(:) = weights./sum(weights);
    
    Nhat = 1/(sum(weights.^2));
    % figure;
    % hold on;
    if Nhat < resample_index
        part_set_new = zeros(size(part_set));
        % perform resampling
        for g=1:particles
            eta = 0.9998*rand(1) + 0.0001;
            % Does not work if it needs every weight to reach eta
            mm1 = 0;
            sumwj = 0;
            while sumwj<=eta
                mm1 = mm1 + 1;
                sumwj = sumwj + weights(mm1);
            end
            m = mm1;
            part_set_new(:, g) = part_set(:, m);
            mhist(k, g) = m;
            % plot(g, m, 'x')
        end
        part_set = part_set_new;
        weights = ones(particles, 1)./particles;
    end
    
        % estimate x and P
    xhat(:, k) = zeros(nx, 1);
    for j=1:particles
        xhat(:, k) = xhat(:, k) + weights(j).*part_set(:, j);
    end
    P(:, :, k) = zeros(nz, nx);
    for j=1:particles
        P(:, :, k) = P(:, :, k) + weights(k)*(part_set(:, j) - xhat(:, k))*(part_set(:, j) - xhat(:, k))';
    end
    
    
end

load('problem4truth.mat');

xtrue = zeros(nx, timesteps);
for t=1:timesteps
    xtrue(:, t) = robot(t).x;
end


if plotting == 1
    figure
    hold on
    plot(1:timesteps, xhat(1, :), 1:timesteps, xtrue(1, :));
    plot(1:timesteps, xhat(1, :) + squeeze(sqrt(P(1, 1, :)))', 'k');
    plot(1:timesteps, xhat(1, :) - squeeze(sqrt(P(1, 1, :)))', 'k');
    title("Comparison of X(1)");
    legend("Xhat(1)", "Xtrue(1)");
    set(gcf, 'position', plotloc1);

    figure
    hold on
    plot(1:timesteps, xhat(2, :), 1:timesteps, xtrue(2, :));
    plot(1:timesteps, xhat(2, :) + squeeze(sqrt(P(2, 2, :)))', 'k');
    plot(1:timesteps, xhat(2, :) - squeeze(sqrt(P(2, 2, :)))', 'k');
    title("Comparison of X(2)");
    legend("Xhat(2)", "Xtrue(2)");
    set(gcf, 'position', plotloc2);
    
    figure
    hold on
    plot(1:timesteps, xhat(3, :), 1:timesteps, xtrue(3, :));
    plot(1:timesteps, xhat(3, :) + squeeze(sqrt(P(3, 3, :)))', 'k');
    plot(1:timesteps, xhat(3, :) - squeeze(sqrt(P(3, 3, :)))', 'k');
    title("Comparison of X(3)");
    legend("Xhat(3)", "Xtrue(3)");
    set(gcf, 'position', plotloc3);

    figure
    hold on
    for i = 2:timesteps
        plot(squeeze(old_parts(i, 1, :)), squeeze(old_parts(i, 2, :)), 'x')
    end
    for i = 2:timesteps
        plot(xtrue(1, i), xtrue(2, i), '.k', 'MarkerSize', 20);
    end
end

set(gcf, 'position', plotloc4);
    xlim([0, 10]);
    ylim([0, 10]);

