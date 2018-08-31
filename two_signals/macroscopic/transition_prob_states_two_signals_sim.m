% Calculates the state transition probabilities for two signals by
% simulations
% pij: for transitions between (i,j) states
clear variables
close all
warning off
set(0, 'defaulttextinterpreter', 'latex');
%%
% Parameters of the system
gz = 15;
N = gz^2;
a0 = 1.5;
rcell = 0.2;
Rcell = rcell*a0;

M_int = [0 1; -1 -1];
Con = [18 16];
Coff = [1 1];
K = [0 12; 13 8];
lambda = [1 1.2];
hill = Inf;
noise = 0;
%}

% Load from file
%{
folder = 'H:\My Documents\Multicellular automaton\app\data\time_evolution\sample_trajectories\typeII';
fname_str = 'two_signal_mult_M_int0_1_-1_-1_period_3_osc-v2';
load(fullfile(folder,fname_str));
N = save_consts_struct.N;
gz = sqrt(N);
a0 = save_consts_struct.a0;
rcell = 0.2;
Rcell = rcell*a0;

M_int = save_consts_struct.M_int;
Con = save_consts_struct.Con;
Coff = save_consts_struct.Coff;
K = save_consts_struct.K;
lambda = [1 save_consts_struct.lambda12];
hill = Inf;
noise = 0;
%}

% simulation parameters
n_trials = 10^3; %10^3; % simulations to consider
delta_t = 1; % # time steps to take

% (1) original hexagonal lattice
%[pos,ex,ey] = init_cellpos_hex(gridsize,gridsize);
%dist = dist_mat(pos,gridsize,gridsize,ex,ey);
% (2) new lattice
mcsteps = 0;
nodisplay = 1;
[pos, dist] = initial_cells_random_markov_periodic(gz, mcsteps, rcell, nodisplay);

% Calculate interaction strength
dist_vec = a0*dist(1,:);
r = dist_vec(dist_vec>0); % exclude self influence
fN = zeros(1, 2);
gN = zeros(1, 2);
for i=1:2
    fN(i) = sum(sinh(Rcell) * sum(exp((Rcell-r)./lambda(i)).*(lambda(i)./r)) ); % calculate signaling strength
    gN(i) = sum(sinh(Rcell)^2 * sum(exp(2*(Rcell-r)./lambda(i)).*(lambda(i)./r).^2 ) ); % calculate noise variance strength
end

% Initial state
iniON = zeros(2);
iniON(1,1) = 50;
iniON(1,2) = 61;
iniON(2,1) = 58;
iniON(2,2) = N - iniON(1,1) - iniON(1,2) - iniON(2,1); %16;

if sum(sum(iniON))~=N
    warning('Wrong probabilities');
end
%% Calculate transition matrix
W = zeros(4, 4); %transition matrix
N_cells_by_state = zeros(1, 4); % count the number of initial cells
t_out = zeros(n_trials, 1);

for trial=1:n_trials
    disp(trial);
    % generate initial lattice

    % subdivide indices for different cell states
    idx = cell(3,1);
    idx_cells_not00 = randperm(N, N-iniON(1,1)); % cell indices of all cells not (0,0)
    part_idx = [1 iniON(1,2) iniON(1,2)+iniON(2,1) iniON(1,2)+iniON(2,1)+iniON(2,2)]; % partition indices
    idx{1} = idx_cells_not00(part_idx(1):part_idx(2));
    idx{2} = idx_cells_not00(part_idx(2)+1:part_idx(3));
    idx{3} = idx_cells_not00(part_idx(3)+1:part_idx(4));

    %{
    %Check that sets have the right number of elems
    check = (numel(idx{1})==iniON(1,2));
    check = check*(numel(idx{2})==iniON(2,1));
    check = check*(numel(idx{3})==iniON(2,2));
    % Check that sets have disjoint elements
    check = check*(numel(union(union(idx{1}, idx{2}), idx{3}))==N-iniON(1,1)); 
    fprintf('Check passed? %d \n', check);
    %}

    % turn genes ON
    cellsIn = zeros(N, 2);
    cellsIn(idx{1}, 2) = 1;
    cellsIn(idx{2}, 1) = 1;
    cellsIn(idx{3}, 1) = 1;
    cellsIn(idx{3}, 2) = 1;
    
    cells = cellsIn;
    cells_hist = {};
    cells_hist{1} = cellsIn;

    % Update cells
    %plothandle = generate_lattice(pos, rcell);
    cell_type = zeros(N, 1);
    disp_mol = 12;
    showI = 1;
    t = 0;
    changed = 1;
    while changed && t < delta_t
        %update_cell_figure_continuum(plothandle, dist, a0, cells, t, disp_mol, showI);
        %waitforbuttonpress;
        %pause(0.5);
        [cellsOut, changed] = update_cells_two_signals_multiply_finite_Hill(cells, dist, M_int, a0,...
                Rcell, Con, Coff, K, lambda, hill, noise);
        cells = cellsOut;   
        cells_hist{end+1} = cells;
        t = t+1;
    end
    t_out(trial) = t;
    
    % Update transition matrix    
    cellsIn_idx = cellsIn*[1; 2];
    cellsOut_idx = cellsOut*[1; 2];
    
    for i=1:N
        W(cellsOut_idx(i)+1, cellsIn_idx(i)+1) = W(cellsOut_idx(i)+1, cellsIn_idx(i)+1) + 1;
    end
    
    % update cell counts
    for i=1:4
        N_cells_by_state(i) = N_cells_by_state(i) + sum(cellsIn_idx==i-1);
    end
    close all
end
W = W./(N_cells_by_state);
% check: all( sum(W, 1) == 1 )

%% Save W
qsave = 0;
folder = 'H:\My Documents\Multicellular automaton\data\two_signals\transition_matrix\sim';
iniON_s = sprintf('%d_%d_%d_%d', iniON(1,1), iniON(1,2), iniON(2,1), iniON(2,2));
a0_s = sprintf('%.2f', a0);
R_s = sprintf('%.2f', rcell);
K_s = sprintf('%d_%d_%d_%d', K(1,1), K(1,2), K(2,1), K(2,2));
Con_s = sprintf('%d_%d', Con(1), Con(2));
M_int_s = sprintf('%d_%d_%d_%d', M_int(1,1), M_int(1,2), M_int(2,1), M_int(2,2));
%I_s = sprintf('%d', mult_dig*round(I, n_dig));
fname_str = strrep(sprintf('W_sim_M_int_%s_N%d_a0_%s_rcell_%s_K_%s_Con_%s_iniON_%s_delta_t_%d_trials_%d',...
    M_int_s, N, a0_s, R_s, K_s, Con_s, iniON_s, delta_t, n_trials), '.', 'p');
if qsave
    save(fullfile(folder, fname_str));
end
%% Plot matrix
h = figure(1);
imagesc(1:4, 1:4, W)
hold on
for i=1:3
    line([0.5, 4.5], [i+0.5, i+0.5], 'Color', 'r', 'LineWidth', 1);
    line([i+0.5, i+0.5], [0.5, 4.5], 'Color', 'r', 'LineWidth', 1);
end
c = colorbar;
set(gca, 'Ydir', 'normal', 'FontSize', 20)
xlabel('$$X(t)$$', 'FontSize', 24)
ylabel('$$X(t+\Delta t)$$', 'FontSize', 24)
c.Label.String = 'Probability';
set(gca, 'XTick', 1:4, 'YTick', 1:4, 'XTickLabels', ...
    {'(0,0)','(1,0)','(0,1)','(1,1)'},...
    'YTickLabels', {'(0,0)','(1,0)','(0,1)','(1,1)'});
set(h, 'Position', [100 100 640 520]);
title(strcat('$$\Delta', sprintf(' t = %d', delta_t), '$$'));
caxis([0 1]);

% save matrix
qsave = 0;
save_figure(h, 10, 8, fullfile(folder, fname_str), '.pdf', qsave);