close all
clear all
warning off

% Calculate the simulated entropy map and save the data

% Parameters of the system
gridsize = 11;
N = gridsize^2;
a0_vec = [0.5 1.5];

% use hexagonal lattice
[pos,ex,ey] = init_cellpos_hex(gridsize,gridsize);
dist = dist_mat(pos,gridsize,gridsize,ex,ey);
dist_vec = dist(1,:);
Son = linspace(1,30,30);
K = linspace(1,20,20);

for idx = 1:numel(a0_vec)
    a0 =a0_vec(idx);
    Rcell = 0.2*a0;
    
    omega_sim = zeros(numel(Son), numel(K));
    omegak = zeros(numel(Son), numel(K), N+1);
    nn = zeros(numel(Son), numel(K), N+1);
    for i=1:numel(Son)
        for j = 1:numel(K)
            [omega_sim(i,j), omegak(i,j,:), nn(i,j,:)] = ...
                sim_entropy(dist, Son(i), K(j), a0, Rcell);
            str_cnt = sprintf('a0: %d --> %.1f%%', idx, ...
                100*(numel(K) * (i - 1) + j)/numel(Son)/numel(K));
            disp(str_cnt)
        end
    end

    fname = strcat(datestr(now,'yyyymmdd_hhMM_'), ...
        sprintf('Son_K_map_N%d_hexagonal.mat', gridsize));
    fname = fullfile(pwd, 'data', fname);
    save(fname);
end
