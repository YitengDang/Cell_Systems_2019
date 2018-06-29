% Simulate pin_pout map for several genetic circuit parameters and save the
% data for future plot. It generates random patterns with order parameter
% I_target

close all
clear variables

data_path = 'D:\eduardopavinat\Dropbox\Matlab codes\data_onecelltype_entropy';

% Parameters of the system
gridsize = 15;
N = gridsize^2;
a0 = 1.5;
Rcell = 0.2*a0;
I_target = 0.2;
save_fig = 0;

% use hexagonal lattice
[dist, pos] = init_dist_hex(gridsize, gridsize);

% points = zeros(10,2);
% points(:,1) = linspace(1,20,10);
% tmp = linspace(1,30,15);
% points(:,2) = tmp(10);

% points(1,:) are the K's and points(2,:) are the Son's
points = zeros(2,2);
points(:,1) = [5, 16];
points(:,2) = [11, 13];

% Runs for all the 
for i = 1:size(points,1)
    disp(i)
    K = points(i,1);
    Son = points(i,2);
    fname = sprintf('pin_poutN%d_Con_%d_K_%d_gz_%d_a0_%d_I_%d.mat', ...
        N, round(Son), round(K), gridsize, 10*a0, round(10*I_target));
    [count, t_av] = count_eq_parallel_I(dist, Son, K, a0, Rcell, I_target);
    % save the data
    fname = fullfile(data_path, 'pin_pout', fname);
    save(fname);
end