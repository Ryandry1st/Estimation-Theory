function [hxg] = littleh(xg,thist)
hxg = zeros(11, 1);
for j=1:11
    hxg(j) = xg(1)*cos(xg(2)*thist(j) + xg(3));
end
end

