%% Initial setup
% Model
% z_j = x1*cost(x2*t_j+x3)+w_j
clear all

nz = 11;
nx = 3;
R = zeros(nz, nz);
for i=1:nz
    for j=1:nz
        if abs(i-j) == 1
            R(i, j) = 0.5;
        elseif i==j
            R(i, j) = 1;
        end
    end
end

thist =[0; 0.1000; 0.2000; 0.3000; 0.4000; 0.5000;
0.6000; 0.7000; 0.8000; 0.9000; 1.0000];

% old zhist
% zhist = [7.7969; 1.4177; -3.0970; -7.6810; -9.8749; -6.1828;
% -0.8212; 4.5074; 8.2259; 9.5369; 6.2827];

% new zhist
zhist = [0.380522980099649
-0.937893887500619
-2.568861909017112
-5.780043794471402
-5.484192626331365
-2.987309363580615
-1.065589854303113
1.286281960892803
3.687753634356679
6.673091483609936
7.148528020590929];

%% Apply the iterations
% start somewhere
xg = [9, 6, 6];
xgbest = xg;
bestcost = P43cost(zhist, thist, xg, R);
nx = 3;
nk = 11;
alpha = 1;
initialcost = P43cost(zhist, thist, xg, R)

    for j=0:20
        for l = 0:20
            xg = [j, 6, l]';
            jold = P43cost(zhist, thist, xg, R);
            for i=1:100
                H = inv(R)*(H_4_3(xg, thist))';
                % How I would have solved it
                % xg = xg + 0.3*(inv(H'*H)*H'*inv(R)*(zhist-xg(1)*cos(xg(2)*thist+xg(3))));
                
                % Use QR factorization
                [Q, R1] = qr(H);
                R1 = R1(1:nx, :);
                Q1 = Q(:, 1:nx);
                lh = littleh(xg, thist);
                delta = (inv(R1)*Q1'*(zhist-lh));
                
                jnew = P43cost(zhist, thist, xg+alpha*delta, R);
                while jnew >= jold +1e-8
                    alpha = alpha/2.0;
                    jnew = P43cost(zhist, thist, xg+alpha*delta, R);
                end
                xg = xg+alpha*delta;
            end
            
            cost = P43cost(zhist, thist, xg, R);
            if cost < bestcost
                xgbest = xg;
                bestcost = cost;
            end
        end
   end

plot(thist, zhist, thist, littleh(xg, thist))
title("Zhist vs Estimator Estimates")
legend("Zhist", "h(xg)")
xlabel("Thist")

H = (H_4_3(xgbest, thist))';
Pxx = inv(H'*inv(R)*H)
xgbest
bestcost

