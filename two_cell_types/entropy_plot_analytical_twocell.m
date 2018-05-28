close all
clear all
warning off

% Calculate the analytical entropy map for a set of parameters ad plot the
% data. It doesn't calculate the Monte Carlo simulation for the entropy.
data_path = '~/Dropbox/Matlab codes/data_twocelltypes_entropy';
% Parameters of the system
gridsize = 11;
N = gridsize^2;
a0 = 1.5;
Rcell = 0.2*a0;
N1 = round(0.6*N);
Nv = [N1 N-N1];

% use hexagonal lattice
[dist, pos] = init_dist_hex(gridsize, gridsize);
dist_vec = dist(1,:);

K1 = 10;
K2 = (1:0.25:30)';
Son1 = 10;
Son2 = 5;

% Calculate analytical entropy
omega = zeros(size(K2));
for i=1:numel(K2)
    omega(i) = entropy_two_cell(dist_vec, Nv, [Son1 Son2], [K1 K2(i)], a0, Rcell);
end

% Plot the result
figure(1)
plot(K2,omega, 'r-')
hold on

K2sim = [1:2:11 20 30];

% Calculate analytical entropy
omegasim = zeros(size(K2sim));
for i=1:numel(K2sim)
    disp(i)
    tic
    omegasim(i) = sim_entropy_parallel_twocell(dist, Nv, [Son1 Son2], [K1 K2sim(i)], a0, Rcell);
    toc
end

plot(K2sim,omegasim, 'ro')
hold off
xlabel('S_{ON}')
ylabel('Entropy')

fname = fullfile(data_path, sprintf('%s_entropy_twocell.mat', ...
    datestr(datetime('now'),'yyyymmddhhMM')));
save(fname);

