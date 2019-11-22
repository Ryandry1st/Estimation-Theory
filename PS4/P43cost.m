function [cost] = P43cost(zhist, thist, xg, R)
cost = sum(abs((zhist-xg(1)*cos(xg(2)*thist+xg(3)))'*inv(R)*(zhist-xg(1)*cos(xg(2)*thist+xg(3))))^(1/2));
end

