close all
clear all
warning off

% Calculate the simulated entropy map and save the data. It runs for
% several grid parameters a0.
% Parameters of the system
gridsize = 15;
N = gridsize^2;
a0_vec = [1.5];
Son = linspace(1,30,30);
K = linspace(1,30,30);

% use hexagonal lattice
[dist, pos] = init_dist_hex(gridsize, gridsize);
dist_vec = dist(1,:);

% Calculate for all a0 in a0_vec
for idx = 1:numel(a0_vec)
    a0 =a0_vec(idx);
    Rcell = 0.2*a0;
    
    % initialize variables
    omega_sim = zeros(numel(Son), numel(K));
    omegak = zeros(numel(Son), numel(K), N+1);
    nn = zeros(numel(Son), numel(K), N+1);
    % Runs for all the combinations of Son and K
    for i=1:numel(Son)
        for j = 1:numel(K)
            [omega_sim(i,j), omegak(i,j,:)] = ...
                sim_entropy_parallel(dist, Son(i), K(j), a0, Rcell);
            str_cnt = sprintf('a0: %.1f --> %.1f%%', a0, ...
                100*(numel(K) * (i - 1) + j)/numel(Son)/numel(K));
            disp(str_cnt) % to monitor the progress
        end
    end
    
    % save the data
    fname = strcat(datestr(now,'yyyymmdd_hhMM_'), ...
        sprintf('Son_K_map_a0%d_N%d_hexagonal.mat', round(10*a0), gridsize));
    fname = fullfile(pwd, 'data','entropy_map', fname);
    save(fname);
end
