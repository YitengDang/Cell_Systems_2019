% Script to calculate the transition matrix in different conditions
% Without I: only takes into account initial p
clear variables
close all
warning off
set(0, 'defaulttextinterpreter', 'latex');

% Parameters
% Hexagonal lattice
gridsize = 11;
[dist, pos] = init_dist_hex(gridsize, gridsize);
N = gridsize^2;

% Feedback
Con = 8;
K = 16;
a0 = 0.5;
rcell = 0.2;
Rcell = 0.2*a0;
M_int = 1; %1: positive, -1: negative

% Load/save folder
fname_str = strrep(sprintf('t_mat_M_int_%d_gz_%d_Con_%.2f_K_%.2f_a0_%.2f_rcell_%.2f', ...
    M_int, gridsize, Con, K, a0, Rcell), '.', 'p');

folder = 'H:\My Documents\Multicellular automaton\data\main\transition_matrix';
fname = fullfile(folder, strcat(fname_str, '.mat'));

if exist(fname_str, 'file') == 2
    load(fname)
else
    t_mat = zeros(N+1);
    for n = 0:N
        disp(n)
        [ptsum, ~, ~] = transition_prob(dist, a0, Rcell, K, Con, n, M_int);
        t_mat(n+1, :) = ptsum;
    end
    save(fname)
end
%%
h=figure(1);
p = (0:N)/N;
imagesc(p, p, t_mat')
c = colorbar;
set(gca, 'Ydir', 'normal', 'FontSize', 20)
xlabel('$$p_{t}$$', 'FontSize', 24)
ylabel('$$p_{t+1}$$', 'FontSize', 24)
c.Label.String = 'Probability';

qsave=1;
if qsave
    folder_fig = 'H:\My Documents\Multicellular automaton\figures\main\transition_matrix';
    fname_fig = fullfile(folder_fig, strcat(fname_str, '_map'));
    save_figure(h, 10, 8, fname_fig, '.pdf')
end