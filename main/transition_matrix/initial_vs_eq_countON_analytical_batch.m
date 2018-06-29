close all
clear variables
warning off

% Runs the branch algorithm for several genetic circuit parameters and save
% the data.

% Parameters of the system
gridsize = 11;
N = gridsize^2;
a0 = 1.5;
Rcell = 0.2*a0;
save_fig = 0;

% Parameters to run
Son_vec = 1:30;
K_vec = 1:20;

% use hexagonal lattice
[pos,ex,ey] = init_cellpos_hex(gridsize,gridsize);
dist = dist_mat(pos,gridsize,gridsize,ex,ey);

for i = 1:numel(Son_vec)
    for j = 1:numel(K_vec)
        Son = Son_vec(i);
        K = K_vec(j);
        % Calculate the transition matrix
        pmatrix = zeros(N+1, N+1);
        pEn = zeros(N+1, 1);
        for n = 0:N
            fprintf('Matrix row: %d\n', n)
            [pmatrix(:,n+1), ~, pEn(n+1)] = transition_prob(dist, a0, Rcell, K, Son, n);
            pmatrix(n+1, n+1) = pmatrix(n+1, n+1) - pEn(n+1); % correct for equilibrium condition
        end
        % Runs the branching algorithm to calculate the map
        out_map = zeros(N+1, N+1);
        t_av = zeros(N+1, 1);
        for n = 0:N
            fprintf('Map row: %d\n', n)
            [out_map(:, n+1), t_av(n+1)] = branching(pmatrix, pEn, 1, n, zeros(N+1, 1), 0, 0);
        end
        fname = fullfile(pwd,'data','pin_pout_analytical',...
            sprintf('grid_%d_Son%d_K%d_a0_%d.mat', gridsize, Son, K, 10*a0));
        save(fname)
    end
end


