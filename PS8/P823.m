%% PS SPF #2
% Unscented Transform
% r = rbar + wr with sigmar
% theta = thetabar + wtheta with sigmatheta
% E[r; theta] = [rbar; thetabar]
% Cov = [sigmar^2, 0; 0, sigmatheta^2]

% transformation is
% x = r*cos(theta)
% y = r*sin(theta)
clear all;
close all;

nx = 2;
nz = 2;

set = 2;

if set==1
    rbar = 76;
    thetabar = -3*pi/180;
    sigmar = 1;
    sigmatheta = pi/180;
else
    rbar = 76;
    thetabar = -3*pi/180;
    sigmar = 1;
    sigmatheta = 15*pi/180;
end

% linearize about rbar and thetabar
xbar = [rbar; thetabar];
% zbar = h(xbar) + 1/2*H*Pxx*H'
dxdr = cos(thetabar);
dxdtheta = -rbar*sin(thetabar);
dydr = sin(thetabar);
dydtheta = rbar*cos(thetabar);

Pxx = [sigmar^2, 0; 0, sigmatheta^2];
H = [dxdr, dxdtheta; dydr, dydtheta];

% zbar = [rbar*cos(thetabar); rbar*sin(thetabar)] + sum(H*Pxx*H', 2);
% truncate after second term
% first term is h(xbar), second is E[Deh] = 0
zbar = [rbar*cos(thetabar); rbar*sin(thetabar)]
Pzz = H*Pxx*H'

%% Unscented Transform
alpha = 10^-3;
beta = 2;
kappa = 0;

chi(:, 1) = xbar;
lambda = alpha^2*(nx+kappa) - nx;
S = chol(Pxx)';

for i=2:nx+1
    chi(:, i) = xbar + sqrt(nx+lambda)*S(:, i-1);
    chi(:, i+nx) = xbar - sqrt(nx+lambda)*S(:, i-1);
end

scriptz = zeros(nz, 2*nx+1);
for j=1:2*nx+1
    rj = chi(1, j);
    thetaj = chi(2, j);
    scriptz(:, j) = [rj*cos(thetaj); rj*sin(thetaj)];
end

zbar2 = lambda/(nx+lambda)*scriptz(:, 1);
for t=2:2*nx+1
    zbar2 = zbar2 + 1/(2*(nx+lambda))*scriptz(:, t);
end

Pzz2 = (lambda/(nx+lambda)+1-alpha^2+beta)*(scriptz(:, 1)-zbar2)*(scriptz(:, 1)-zbar2)';
for t=2:2*nx+1
    Pzz2 = Pzz2 + 1/(2*(nx+lambda))*(scriptz(:, t)-zbar2)*(scriptz(:, t)-zbar2)';
end

zbar2
Pzz2
% diff_zbar = zbar-zbar2
% diff_Pzz = Pzz-Pzz2

%% Generate vector and compare
ton = 1e6;

tonofXs = GaussGen(xbar, Pxx, nx, ton);
% % samp_cov = zeros(nx, nx);
% % for j=1:ton
% %     samp_cov = samp_cov + tonofXs(:, j)*tonofXs(:, j)';
% % end
% % 
% % samp_cov = samp_cov./ton - xbar*xbar'

tonofZs = zeros(size(tonofXs));

for j=1:ton
    rj = tonofXs(1, j);
    thetaj = tonofXs(2, j);
    tonofZs(:, j) = [rj*cos(thetaj); rj*sin(thetaj)];
end

samp_mean = zeros(nx, 1);
for j=1:ton
    samp_mean = samp_mean + tonofZs(:, j);
end
samp_mean = samp_mean ./ ton

samp_cov = zeros(nx, nx);
for j=1:ton
    samp_cov = samp_cov + tonofZs(:, j)*tonofZs(:, j)';
end

samp_cov = samp_cov./ton;
samp_cov = samp_cov - samp_mean*samp_mean'

% parta_zbar_diff = zbar - samp_mean
% partb_zbar_diff = zbar2 - samp_mean
% 
% % Pxx does not seem to change too much..??
% parta_Pzz_diff = Pzz - samp_cov
% partb_Pzz_diff = Pzz2 - samp_cov

