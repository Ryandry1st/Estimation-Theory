% Solve P4-6
% note x = 0 and z = 0 are the truth
% start with xhat = 1.5
% h(x) = atan(x)
% J = ||z-h(x)||
clear all;
stepper = 0;

alpha = 0.92;
runs = 50;
xg = 1.5;
z = 0;
oldxg = zeros(runs, 3);

alphas = 0:1/1e6:1;

for j=1:3
    for i=1:runs
        h = atan2(xg, 1);

        J0 = abs(norm(0-h));
        H = cos(h)^2;
        dxg = alpha * (H*H)^-1*H*(-h);
        % No step size adjustment
        if stepper == 0
            oldxg(i, 1) = xg;
            xg = xg + dxg;
        elseif stepper == 1
            J1 = abs(norm(0 - atan(xg+dxg)));
            while J1 >= J0
                alpha = alpha/2.0;
                dxg = alpha * (H*H)^-1*H*(-h);
                J1 = abs(norm(0 - atan(xg+dxg)));
            end
            oldxg(i, 2) = xg;
            xg = xg + dxg;
        else
            dxgs = alphas.*inv((H*H))*H*(-h);
            Js = abs(atan(xg+dxgs));
            [argval, argmin] = min(Js);
            oldxg(i, 3) = xg;
            xg = xg + dxgs(argmin);
        end
        % xg

    end
    if j==1
        stepper=1;
    else
        stepper=2;
    end
    xg = 1.5;
    alpha = 1.0;
end

x = -6:1/100:6;

figure
plot(x, abs(atan2(x, 1)))
hold on
plot(oldxg(:, 1), abs(atan2(oldxg(:, 1), 1)))
plot(oldxg(:, 2), abs(atan(oldxg(:, 2))))
plot(oldxg(:, 3), abs(atan(oldxg(:, 3))))
legend("Cost Function", "Alpha=0.92", "Alpha step adj", "My method")
title("Comparison of step size adjustment for GN")

min(oldxg(:, 1))
min(oldxg(:, 2))
min(oldxg(:, 3))

