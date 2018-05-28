% Script to calculate the transition matrix in different conditions
clear variables
close all
warning off

%data_path = 'C:\Users\eduardopavinat\Dropbox\Matlab codes\data_onecelltype_entropy';
folder = 'transition_matrix';

% Parameters
% Hexagonal lattice
gridsize = 11;
[dist, pos] = init_dist_hex(gridsize, gridsize);
N = gridsize^2;
% Feedback
Con = 15;
K = 6;
a0 = 1.5;
Rcell = 0.2*a0;
% Spatial order
I = 0.5;

filename = sprintf('t_mat_grid%d_Con%d_K%d_a0%d_R%d_I%d.mat', ...
    gridsize, round(Con), round(K), round(10*a0), ...
    round(100*Rcell), round(100*I));

fname = fullfile(pwd, 'data', folder, filename);

if exist(filename, 'file') == 2
    load(fname)
else
    t_mat = zeros(N+1);
    for n = 0:N
        disp(n)
        [ptsum, ~, ~] = transition_prob_with_I(dist, a0, Rcell, K, Con, n, I);
        t_mat(n+1, :) = ptsum;
    end
    save(fname)
end

p = (0:N)/N;
imagesc(p, p, t_mat')
c = colorbar;
set(gca, 'Ydir', 'normal')
xlabel('p_{t}', 'FontSize', 24)
ylabel('p_{t+1}', 'FontSize', 24)
