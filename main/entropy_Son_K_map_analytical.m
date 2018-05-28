close all
clear all
warning off

% Calculate the analytical entropy map and save the data. It runs for
% several a0 specified in a0_vec.

% Parameters of the system
gridsize = 30;
N = gridsize^2;
a0_vec = [0.5 1.5];
Son = linspace(1,30,30);
K = linspace(1,20,20);

% use hexagonal lattice
[dist, pos] = init_dist_hex(gridsize, gridsize);
dist_vec = dist(1,:);

% Calculate for all a0 in a0_vec
for idx = 1:numel(a0_vec)
    a0 =a0_vec(idx);
    Rcell = 0.2*a0;
    % initialize variables
    omega_sim = zeros(numel(Son), numel(K));
    % Runs for all the combinations of Son and K
    for i=1:numel(Son)
        for j = 1:numel(K)
            omega_sim(i,j) = ...
                entropy_eq_sphere(dist_vec, Son(i), K(j), a0, Rcell);
            str_cnt = sprintf('a0: %d --> %.1f%%', idx, ...
                100*(numel(K) * (i - 1) + j)/numel(Son)/numel(K));
            disp(str_cnt) % to monitor the progress
        end
    end

    % save the data
    fname = strcat(datestr(now,'yyyymmdd_hhMM_'), ...
        sprintf('Son_K_map_N%d_a0%dhexagonal_analytical.mat', gridsize, 10*a0));
    fname = fullfile(pwd, 'data', fname);
    save(fname);
end
