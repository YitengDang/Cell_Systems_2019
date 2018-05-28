%% Calculates I for a uniform lattice and plots it as a function of p
syms x
p1 = x^2 - (2*x-1)^2;
p2 = 4*p*(1-p);
[q, r] = quorem(p1, p2, x);
%%
fplot(q, [0 1])