clear variables
close all

gridsize = 15;
a0 = 0.5;
Rcell = 0.2*a0;

% initialize cells in grid
[dist, pos] = init_dist_hex(gridsize, gridsize);

% Matrix used to calculate the influence of cells on each other
M = zeros(size(dist));
M(dist > 0) = exp(Rcell - a0*dist(dist > 0))*sinh(Rcell)./(a0*dist(dist > 0));

cov_Y = M'*M;
gN = cov_Y(1,1);

nsamples = 10000;
out = zeros(1, nsamples);
N = size(M,1);
for i = 1:nsamples
    out(i) = sum(mvnrnd(zeros(1, N), cov_Y))/N;
end

mean(out)
std(out)
histogram(out)