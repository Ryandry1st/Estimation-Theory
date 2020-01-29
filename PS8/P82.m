%% PS SPF # 1
clear all;
close all;

proof = 0;

if proof == 1
    nx = 1000;
    z = randn(nx, 1);

    % x = A*z + b ~ N(b, AA')
    % Proof
    A = rand(nx, nx);
    b = ones(nx, 1);
    x = A*z + b;
    mean_diff = mean(x - b);
    cov_diff = norm(mean((x-b)*(x-b)' - A*z*z'*A'));

    % Use eigenvalue decomposition to solve for A for Pxx
    Pxx = rand(nx, nx);
    xbar = b;
    [V, D] = eig(Pxx);
    A = V*sqrt(D);
else
    nx = 3;
    n = 1000000;
    xbar = ones(nx, 1);
    Pxx = eye(nx) + abs(rand(nx, nx));
    
    xbig = GaussGen(xbar, Pxx, nx, n);
    xbig_mean = mean(xbig, 2);
    avg_mean_diff = xbig_mean - xbar
    for k=1:n
        cov_diff(:, :, k) = (xbig(:, k)-xbar)*(xbig(:, k)-xbar)'-Pxx;
    end
    avg_cov_diff = norm(mean(cov_diff, 3))
    
end



