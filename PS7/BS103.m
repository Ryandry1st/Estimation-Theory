%% BS 10.3
% rho = sqrt(x1^2 + x2^2)
% theta = arctan(x2/x1)
% z = [z1; z2] = [rho+w1; theta + w2] = h(x1, x2) + w
% R = diag(sigmarho^2, sigmatheta^2)
clear all

% y = inv(h(z) = [x1; x2] + Wc
% y = [z1*cos(z2); z1*sin(z2)]

% Find Rc using linearization at rho=10^5, theta=45
rbar = 10^5;
thetabar = 45/180*pi; % convert to radians

% linearize 
zbar = [rbar; thetabar];
dxdr = cos(thetabar);
dxdtheta = -rbar*sin(thetabar);
dydr = sin(thetabar);
dydtheta = rbar*cos(thetabar);

sigmap = 100;
sigmatheta = 0.5/180*pi;
Pzz = [sigmap^2, 0; 0, sigmatheta^2];
H = [dxdr, dxdtheta; dydr, dydtheta];

ybar = [rbar*cos(thetabar); rbar*sin(thetabar)]

Pxx = H*Pzz*H'

N = 10000;
nx = 2;
R = diag([100^2, (0.5/180*pi)^2]);
w = GaussGen(0, R, nx, N);
z = zbar + w;
y = zeros(nx, N);
for i=1:N
    y(:, i) = [z(1, i)*cos(z(2, i)); z(1, i)*sin(z(2, i))];
end
y_avg = mean(y, 2)
Wcmean = mean(y-ybar, 2)


cov_test = zeros(nx, nx);
cov_test = cov(y')