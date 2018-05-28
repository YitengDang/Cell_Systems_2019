close all
clear variables

% Simulate pin_pout map for several genetic circuit parameters and save the
% data for future plot.

% Path to save the data
data_path = 'C:\Users\eduardopavinat\Dropbox\Matlab codes\data_onecelltype_entropy';

% Parameters of the system
gridsize = 15;
N = gridsize^2;
a0 = 1.5;
Rcell = 0.2*a0;
save_fig = 0;

% points(1,:) are the K's to be tested and points(2,:) are the Son's to be
% tested

points = [6 15; 12 15; 20 15];

% points = zeros(10,2);
% points(:,1) = linspace(1,20,10);
% tmp = linspace(1,30,15);
% points(:,2) = tmp(10);

% use hexagonal lattice
[pos,ex,ey] = init_cellpos_hex(gridsize,gridsize);
dist = dist_mat(pos,gridsize,gridsize,ex,ey);

% Runs for all points
for i = 1:size(points,1)
    disp(i)
    K = points(i,1);
    Son = points(i,2);
    fname = sprintf('pin_poutN%d_Con_%d_K_%d_gz_%d_a0_%d.mat', ...
        N, round(Son), round(K), gridsize, 10*a0);
    [count, t_av, I_av] = count_eq_parallel(dist, Son, K, a0, Rcell);
    fname = fullfile(data_path, 'pin_pout', fname);
    save(fname);
end