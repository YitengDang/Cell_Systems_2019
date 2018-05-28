close all
clear variables

% Simulate pin_pout map for several genetic circuit parameters and save the
% data for future plot

data_path = 'D:\eduardopavinat\Dropbox\Matlab codes\data_onecelltype_entropy';

% Parameters of the system
gridsize = 15;
N = gridsize^2;
a0 = 1.5;
Rcell = 0.2*a0;
save_fig = 0;

% use hexagonal lattice
[dist, pos] = init_dist_hex(gridsize, gridsize);

for b = 4:4
    switch b
           
        case 1
            points = zeros(15,2);
            points(:,2) = linspace(1,30,15);
            tmp = linspace(1,20,10);
            points(:,1) = tmp(3);
            alpha = 0.05*points(:,1);

        case 2
            points = zeros(15,2);
            points(:,2) = linspace(1,30,15);
            tmp = linspace(1,20,10);
            points(:,1) = tmp(8);
            alpha = 0.05*points(:,1);
            
        case 3
            points = zeros(10,2);
            points(:,1) = linspace(1,20,10);
            tmp = linspace(1,30,15);
            points(:,2) = tmp(5);
            alpha = 0.05*points(:,1);
            
        case 4
            points = zeros(10,2);
            points(:,1) = linspace(1,20,10);
            tmp = linspace(1,30,15);
            points(:,2) = tmp(10);
            alpha = 0.05*points(:,1);
            
            
    end
    for i = 1:size(points,1)
        disp(i)
        K = points(i,1);
        Son = points(i,2);
        noise = alpha(i);
        fname = sprintf('pin_poutN%d_Con_%d_K_%d_gz_%d_a0_%d_noise_%d.mat', ...
            N, round(Son), round(K), gridsize, round(10*a0), round(10*noise));
        [count, t_av, I_av] = count_eq_parallel_noise(dist, Son, K, a0, Rcell, noise);
        fname = fullfile(data_path, 'pin_pout_noise', fname);
        save(fname);
    end
end