function [cost, x, y, ra_hat, rb_hat, ta_hat, tb_hat] = Cost4_4(xg, thist, ra, rb, ta, tb, nk)

Rs = diag([10^2*ones(1, nk), 30^2*ones(1, nk)]);

tb_hat = 0;
ta_hat = 0;

x = xg(1) + xg(3).*thist;
y = xg(2) + xg(4).*thist - 4.9*thist.^2;

ra_hat = sqrt((4.1e5 - x).^2 + y.^2);
rb_hat = sqrt((4.4e5 - x).^2 + y.^2);

if isscalar(ta)
    r_hat = [ra_hat; rb_hat];
    r = [ra; rb];
else
    ta_hat = atan2(y, -x + 4.1e5);
    tb_hat = atan2(y, -x + 4.4e5);
    
    Rs = diag([10^2*ones(1, nk), 30^2*ones(1, nk), 0.01^2*ones(1, nk), 0.03^2*ones(1, nk)]);
    r_hat = [ra_hat; rb_hat; ta_hat; tb_hat];
    r = [ra; rb; ta; tb];
end



cost = sqrt((r-r_hat)'*inv(Rs)*(r-r_hat));
end

