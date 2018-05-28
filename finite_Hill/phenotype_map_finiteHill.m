% Calculate the simulated phenotype score map with finite Hill coefficient and save the data. It runs for
% several grid parameters a0. Saves data automatically.
% Adapted from entropy_Son_K_map.m
close all
clear all
warning off
%%
% Parameters of the system
gridsize = 15;
N = gridsize^2;
a0_vec = [1.5];
Con = linspace(2,30,50);
K = linspace(1,20,50);
hill = 1;
initialID = 'uniform'; %ID for initial configuration
%Examples: 'uniform', 'allON', 'singleOFF'

% use hexagonal lattice
[dist, pos] = init_dist_hex(gridsize, gridsize);
dist_vec = dist(1,:);

% Calculate for all a0 in a0_vec
for idx = 1:numel(a0_vec)
    a0 =a0_vec(idx);
    Rcell = 0.2*a0;
    
    % initialize variables
    xmean = zeros(numel(Con), numel(K));
    xstd = zeros(numel(Con), numel(K));
    tmean = zeros(numel(Con), numel(K));

    % Runs for all the combinations of Son and K
    for i=1:numel(Con)
        for j = 1:numel(K)
            [xmean(i,j), xstd(i,j), tmean(i,j)] = phenotype_sim(dist,...
                Con(i), K(j), a0, Rcell, hill);
            str_cnt = sprintf('a0: %.1f --> %.1f%%', a0, ...
                100*(numel(K) * (i - 1) + j)/numel(Con)/numel(K));
            disp(str_cnt) % to monitor the progress
        end
    end
    
    % save the data
    i=1;
    fname_str = strrep(sprintf('Con_K_map_a0%d_N%d_hill%.1f_hexagonal_',...
        round(10*a0), gridsize^2, hill), '.', 'p');
    fname = fullfile(pwd, 'data','phenotype_finite_Hill',...
            strcat(fname_str, initialID,'-v',int2str(i),'.mat'));
    while exist(fname, 'file')==2
        i=i+1;
        fname = fullfile(pwd, 'data','phenotype_finite_Hill',...
            strcat(fname_str, initialID, '-v',int2str(i),'.mat'));
    end
    save(fname);
    
end