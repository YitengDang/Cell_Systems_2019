% Time evolution of a system with visualization of the dynamics without
% noise showing the count of nearest neighbors that are ON
close all
clear all
warning off
set(0, 'defaulttextinterpreter', 'latex');

%%
% lattice parameters
gz = 16;
N = gz^2;
a0 = 1.5;
rcell = 0.2;
Rcell = rcell*a0;
lambda12 = 1.2;
lambda = [1 lambda12];

% circuit parameters
M_int = [0 1; -1 1];
Con = [18 16];
Coff = [1 1];
K = [0 9; 11 4];
hill = Inf;
noise = 0;

% K for all cells
K_all = zeros(N, 2, 2);
for i=1:N
    K_all(i, :, :) = K;
end

% initial conditions
p0 = [0.5 0.5];
iniON = round(p0*N);
%iniON = 1

% Initialize parameters
%[pos,ex,ey] = init_cellpos_hex(gridsize,gridsize);
%dist = dist_mat(pos,gridsize,gridsize,ex,ey);
mcsteps = 0;
[pos, dist, ~, ~] = initial_cells_random_markov_periodic(gz, mcsteps, rcell);

dist_vec = a0*dist(1,:);
r = dist_vec(dist_vec>0); % exclude self influence
fN = sum(sinh(Rcell)*sum(exp(Rcell-r)./r)); % calculate signaling strength

% generate cell_type (0 case type 1, 1 case type 2)
cell_type = zeros(N,1); % all the same here

%% Plot phase diagram
%M_int = 1; lambda12 = 1;
h = plot_phase_diagram(gz, a0, rcell, M_int, K, Con, lambda12);

%% K profile
% Choose which interaction to make spatially dependent
% interaction = [i j] means the j->i interaction
int_wave = [2 1];

% Parameters
Lx = 1; %default size

Ax = 0.6;
nx = 1; %number of times the wave fits into x-range
lambda_x = 1/nx*Lx;
Ay = 0.0;
ny = 1;
lambda_y = 1/ny*sqrt(3)/2*Lx;
wave_x = Ax.*square(pos(:, 1).*(2*pi/lambda_x));
wave_y = Ay.*square(pos(:, 2).*(2*pi/lambda_y));

% sinusoidal wave
% K_func = @(x, y) K.*(1 + Ax.*sin(2.*pi.*x/lambda_x) + Ay.*sin(2.*pi.*y/lambda_y));
% K_all = K_func(pos(:, 1), pos(:, 2));

% square wave
K_all(:, int_wave(1), int_wave(2)) = K(int_wave(1),int_wave(2)).*(1 + wave_x + wave_y);

%{
figure;
x = 0:0.01:1;
y = 0:0.01:1;
plot(pos(:, 1), K_all, 'bo');
%}

%{
figure;
imagesc(reshape(K_cells_all, gz, gz))
colorbar;
set(gca, 'YDir', 'normal');
%}
%% Check that K values should still be able to generate waves (nearest-neighbour calculation)
Kmin = min(K_all(:, int_wave(1), int_wave(2)));
Kmax = max(K_all(:, int_wave(1), int_wave(2)));

Kmin_set = K; 
Kmin_set(int_wave(1), int_wave(2)) = Kmin;
Kmax_set = K; 
Kmax_set(int_wave(1), int_wave(2)) = Kmax;

conditions_met = zeros(2, 3);
wave_types = {'plane', 'inward', 'outward'};
for i=1:numel(wave_types)
    wave_type = wave_types{i};
    conditions_met(1, i) = all(trav_wave_conditions_nn(dist, Con, Kmin_set,...
        rcell, a0, lambda, wave_type));
    conditions_met(2, i) = all(trav_wave_conditions_nn(dist, Con, Kmax_set,...
        rcell, a0, lambda, wave_type));
end

disp('Nearest neighbour conditions for plane waves met?');
disp('Rows: [Kmin; K_max], Columns = [plane, inward, outward]');
disp(conditions_met);

%% Plot profile
h = figure;
% all sizes in units of pixels
Sx = 800; %ax.Position(3); %512;
Sy = (sqrt(3)/2*(gz-1)+2)/(3/2*(gz+1))*Sx;
a0_fig = Sx/(sqrt(N)+1); %Lx/(3/2*(gz+1));
Rcell_fig = rcell*a0_fig;
set(h, 'Position', [100 100 Sx Sy]);

% set image properties
set(gca, 'YTick', [], 'XTick', [], 'Color', [0.8 0.8 0.8]);
title(gca, sprintf('$$A_x = %.1f, A_y = %.1f$$', Ax, Ay), 'FontSize', 20);
Lx = 1;
d = 2*rcell*Lx/(sqrt(N)+1);
Ly = sqrt(3)/2*Lx;
xlim([-d Lx+d]);
ylim([-d Ly+d]);

% --plot cells--
hold on
% colours
c_all = K_all(:, int_wave(1), int_wave(2));
clr_k = zeros(N, 3); % black boundaries
%markers = {'o', 's'};

scatter(pos(:,1), pos(:,2), Rcell_fig^2, c_all, 'filled', 'o');
scatter(pos(:,1), pos(:,2), Rcell_fig^2, clr_k, 'o'); % plot cell boundaries

% Plot box outline
plot([0 Lx], [0 0], 'k--');
plot([0 Lx], [Ly Ly], 'k--');
plot([0 0], [0 Ly], 'k--');
plot([Lx Lx], [0 Ly], 'k--');

% colorbar
c = colorbar;
%c.Ticks = 0:0.2:1;
set(c, 'FontSize', 12);
c.Label.String = '$$K(i)$$';
ylabel(c, sprintf('$$K^{(%d %d)}$$', int_wave(1), int_wave(2)),...
    'Interpreter', 'latex', 'FontSize', 16)
map = 'parula';
colormap(map);
%pause(0.5);

% Save map
qsave = 0;
colored_background = 1;
folder = 'H:\My Documents\Multicellular automaton\figures\parameter_gradient';
fname_str = strrep(sprintf('gradient_profile_cells_Ax_%.1f_Ay_%.1f_nx_%d_ny_%d',...
    Ax, Ay, nx, ny), '.','p');
save_figure(h, 0, 0, fullfile(folder, fname_str), '.pdf', qsave, colored_background)

%%
% initialize ON cells
cells = zeros(N,2);
cells(randperm(N,iniON(1)), 1) = 1;
cells(randperm(N,iniON(2)), 2) = 1;

% variables
t = 0;
I = [];
Non = [];
h = [];
cells_hist = {};

% initialize figure
hin = figure;
disp_mol = 12;
showI = 0;
plot_handle = reset_cell_figure(hin, pos, rcell);
update_figure_periodic_scatter(plot_handle, cells, t, disp_mol, showI, a0, dist)

% save vars and update cells
cells_hist{end+1} = cells;
[cells_out, changed] = ...
    update_cells_two_signals_multiply_v2(cells, dist, M_int, a0, Rcell, Con, Coff, K_all, lambda, noise);
while changed
    pause(0.1);
    t = t+1;
    update_figure_periodic_scatter(plot_handle, cells_out, t, disp_mol, showI, a0, dist)
    cells_hist{end+1} = cells_out;
    cells = cells_out;
    [cells_out, changed] = ...
    update_cells_two_signals_multiply_v2(cells, dist, M_int, a0, Rcell, Con, Coff, K_all, lambda, noise);
end

%% Save result
%{
data_path = 'H:\My Documents\Multicellular automaton\data\two_signals\parameter_gradient';

% variables
positions = pos;
distances = dist;
t_out = numel(cells_hist)-1;
I0 = [Inf Inf];
sim_ID = 'two_signal_mult';
hill_ID = 'hill_Inf';
ini_ID = '';
I_ini_str = '';
save_vars = {N, a0, M_int, K, Con, hill, noise, p0, I0, rcell, lambda12...
    sim_ID, hill_ID, ini_ID, I_ini_str, mcsteps};
save_vars_lbl = {'N', 'a0', 'M_int', 'K', 'Con', 'hill', 'noise', 'p_ini', 'I_ini',...
    'rcell', 'lambda12', 'sim_ID', 'hill_ID', 'ini_ID', 'I_ini_str', 'mcsteps'};
save_consts_struct = cell2struct(save_vars, save_vars_lbl, 2);

save_gradient_vars = {int_wave, Ax, Ay, lambda_x, lambda_y, Lx};
save_gradient_vars_lbl = {'int_wave', 'Ax', 'Ay', 'lambda_x', 'lambda_y', 'Lx'};
save_gradient_struct = cell2struct(save_gradient_vars, save_gradient_vars_lbl, 2);

fname_str = sprintf('Parameter_gradient_K_%d_%d_square_wave', ...
    int_wave(1), int_wave(2));
i = 1;
fname = fullfile(data_path, ...
    strcat(fname_str,'-v',int2str(i),'.mat'));
while exist(fname, 'file') == 2
    i=i+1;
  fname = fullfile(data_path, ...
      strcat(fname_str,'-v',int2str(i),'.mat'));
end
    
save(fname, 'save_consts_struct',...
    'cells_hist', 't_out', 'changed', 'positions', 'distances', 'save_gradient_struct');
%}
%%
%{
data_path = 'H:\My Documents\Multicellular automaton\data\two_signals\parameter_gradient';
load(fullfile(data_path, 'Parameter_gradient_K_2_1_square_wave-v1'))
%}