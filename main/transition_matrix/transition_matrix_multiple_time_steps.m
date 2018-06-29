% Scrpt to calculate the transition matrix in different conditions
% Without I: only takes into account initial p
% Multiple time_steps: evolve the systme over multiple time steps
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
Rcell = rcell*a0;
M_int = 1; %1: positive, -1: negative

% Time steps to advance
tsteps = 1;

% Load/save folder
fname_str = strrep(sprintf('t_mat_M_int_%d_gz_%d_Con_%.2f_K_%.2f_a0_%.2f_rcell_%.2f_tsteps_%d', ...
    M_int, gridsize, Con, K, a0, rcell, tsteps), '.', 'p');

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
    % repeat matrix multiplication
    t_mat_mult = mpower(t_mat, tsteps);
    save(fname)
end

%%
h=figure(1);
p = (0:N)/N;
imagesc(p, p, t_mat_mult')
c = colorbar;
set(gca, 'Ydir', 'normal', 'FontSize', 20)
xlabel('$$p_{t}$$', 'FontSize', 24)
ylabel(sprintf('$$p_{t+%d}$$', tsteps), 'FontSize', 24)
c.Label.String = 'Probability';

qsave=1;
if qsave
    folder_fig = 'H:\My Documents\Multicellular automaton\figures\main\transition_matrix';
    fname_fig = fullfile(folder_fig, strcat(fname_str, '_map'));
    save_figure(h, 10, 8, fname_fig, '.pdf')
end

%% 
h2=figure(2);
plot(p, diag(t_mat_mult), 'o-');
xlabel('$$p_{t}$$');
ylabel(sprintf('$$P(p_{t+%d} = p_{t})$$', tsteps));
set(gca,'FontSize', 20)
xlim([0 1]);
ylim([0 1]);

qsave=1;
if qsave
    folder_fig = 'H:\My Documents\Multicellular automaton\figures\main\transition_matrix';
    fname_fig = fullfile(folder_fig, strcat(fname_str, '_P_equal') );
    save_figure(h2, 10, 8, fname_fig, '.pdf')
end
