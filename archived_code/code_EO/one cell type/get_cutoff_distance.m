% Plot cutoff distance. This is calculated in two different ways. First we
% assume that at a certain distance the sum of all the remaining cells (ON)
% positioned at that distance cannot give the contribution of Soff = 1. In
% the second we assume that the sum of all remaining cells (ON) positioned
% at the cutoff distance cannot promote transition (they cannot reach K-1).

% I am not sure if this approach is valid (Eduardo 27-7-16).
clear variables

% Parameters
Son_vec = 2:0.2:30;
K_vec = 2:0.2:20;

a0 = 1;
Rcell = 0.2*a0;
N = 121;

out = zeros(numel(Son_vec), numel(K_vec));
out2 = zeros(numel(Son_vec),1);

for i = 1:numel(Son_vec)
    for j = 1:numel(K_vec)
        factor = (K_vec(j)-1)/(N-1)/Son_vec(i)/Rcell/exp(Rcell);
        f = @(x) exp(-x) - factor*x;
        out(i,j) = fzero(f, a0);
    end
    factor = 1/(N-1)/Son_vec(i)/Rcell/exp(Rcell);
    f = @(x) exp(-x) - factor*x;
    out2(i) = fzero(f, a0);
end