%% Solve for Gk

syms tk tkp1 a b x dt
A = [0, 1; 0, -a];
B = [0; b];

Gk = int(expm(A*(tk+dt - x))*B, x, [tk, tk+dt])