% Addition to PS4 #1
% generate all possible initial values
clear all;
close all;

x1 = -10:10;
x2 = -10:10;

% total possible convergences
final = zeros(21, 21, 2);
uniquevalues = [];

% run each possible input and determine if it converges
for i=-10:10
    for j=-10:10
        xg = [i; j];
        % ensure convergence
        for k=1:20
            xg = P4_1_grad_desc(xg);
        end
        cost = norm([xg(1) + xg(2) + xg(1)*xg(2)+5; xg(1)^2+2*xg(2)-xg(2)^2-2], 2);
        if cost < 1e-5
            final(i+11, j+11, :) = round(xg, 4); % for determining uniqueness use rounding
            if ~ ismember(final(i+11, j+11, :), uniquevalues)
                uniquevalues = [uniquevalues, final(i+11, j+11, :)];
            end
        end
        % note that [0; 0] is not a convergent spot so any of these are not
        % convergent
        
    end
end

% determine uniqueness
% plot which values are convergent
uniquevalues = squeeze(uniquevalues);

figure
plot(uniquevalues(1, 1), uniquevalues(1, 2) , 'rx', 'MarkerSize',12);
hold on
title("Initial Guesses and Convergent Points")
xlim([-10, 10])
ylim([-10, 10])
xlabel("X1")
ylabel("X2")

plot(uniquevalues(2, 1), uniquevalues(2, 2), 'ro', 'MarkerSize',12);
hold on

nonconv = 0;

for i=-10:10
    for j=-10:10
        if squeeze(final(i+11, j+11, :)) == uniquevalues(1, :)'
            plot(i, j, 'bx')
            hold on
        elseif squeeze(final(i+11, j+11, :)) == uniquevalues(2, :)'
            plot(i, j, 'bo')
            hold on
        else
            nonconv = nonconv + 1;
        end
    end
end

nonconv