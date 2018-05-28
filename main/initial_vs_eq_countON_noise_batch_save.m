close all
clear variables

% Simulate pin_pout map for several genetic circuit parameters and save the
% data for future plot. This uses the exact simulation to generate the map.
% This uses the exact simulation to generate the map with noise in the 
% dynamics

% Parameters of the system
gridsize = 15;
N = gridsize^2;
%a0 = 1.5;
%Rcell = 0.2*a0;
save_fig = 0;

% use hexagonal lattice
[pos,ex,ey] = init_cellpos_hex(gridsize,gridsize);
dist = dist_mat(pos,gridsize,gridsize,ex,ey);

% points(1,:) are the K's to be tested and points(2,:) are the Son's to be
% tested

% points = zeros(15,2);
% points(:,2) = linspace(1,30,15);
% tmp = linspace(1,20,10);
% points(:,1) = tmp(8);
% alpha = 0.1*points(:,1); % noise vector
%K = 6;
%Con = 15;
alpha = [0.8 1];
points = [8 15 0.5; 20 15 1.5];
%points = repmat([K Con], numel(alpha), 1);

%%
for i = 1:size(points,1)
    disp(i)
    K = points(i,1);
    Con = points(i,2);
    a0 = points(i,3);
    Rcell = a0*0.2;
    noise = alpha(i);

    [count, t_av, I_av] = count_eq_parallel_noise(dist, Con, K, a0, Rcell, noise);
    
    % Save the data
    fname_str = sprintf('pin_pout_N%d_Con_%d_K_%d_gz_%d_a0_%d_noise_%d.mat', ...
        N, round(Con), round(K), gridsize, round(10*a0), round(10*noise));
    fname = fullfile(pwd, 'data', 'pin_pout', 'noise', fname_str);
    save(fname);
end