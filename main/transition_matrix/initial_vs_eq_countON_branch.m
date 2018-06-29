close all
clear variables
warning off

%fdir = 'C:\Users\eduardopavinat\Dropbox\Matlab codes\data_onecelltype_entropy\pin_pout_branch';

% Calculate the p_in p_eq map using the branch algorithm for a single set 
% of parameters

% Parameters of the system
gridsize = 15;
N = gridsize^2;
a0 = 0.5;
Rcell = 0.2*a0;
save_fig = 1;
Son = 14;
K = 19;

% use hexagonal lattice
[pos,ex,ey] = init_cellpos_hex(gridsize,gridsize);
dist = dist_mat(pos,gridsize,gridsize,ex,ey);

% Calculate the transition matrix for the branching algorithm
pmatrix = zeros(N+1, N+1);
pEn = zeros(N+1, 1);
for n = 0:N
    fprintf('Matrix row: %d\n', n)
    [pmatrix(:,n+1), ~, pEn(n+1)] = transition_prob(dist, a0, Rcell, K, Son, n);
    pmatrix(n+1, n+1) = pmatrix(n+1, n+1) - pEn(n+1); % correct for equilibrium condition
end
% Calculate the map
out_map = zeros(N+1, N+1);
t_av = zeros(N+1, 1);
for n = 0:N
    fprintf('Map row: %d\n', n)
    [out_map(:, n+1), t_av(n+1)] = branching(pmatrix, pEn, 1, n, zeros(N+1, 1), 0, 0);
end

% Plot the map
p = (0:N)/N;

figure(1)
imfig = imagesc(p,p,out_map);
set(imfig, 'AlphaData', out_map>1e-10)
set(gca, 'ydir', 'normal')

% Plot th average number of steps
figure(2)
plot(p, t_av)

% Save the data
fname = fullfile(pwd,'dynamics','pin_pout_branch',...
    sprintf('grid_%d_Son%d_K%d_a0_%d.mat', gridsize, Son, K, 10*a0));
save(fname)

