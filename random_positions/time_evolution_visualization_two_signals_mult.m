close all
clear all
% maxNumCompThreads(4);
% warning off
%% (1) input parameters
%{
% lattice parameters
gz = 15;
N = gz^2;
a0 = 1.5;
rcell = 0.2;
Rcell = rcell*a0;

% circuit parameters
Con = [18 16];
Coff = [1 1];
M_int = [1 1; -1 -1];
K = [3 12; 13 20]; % K(i,j): sensitivity of type i to type j molecules
lambda = [1 1.2]; % diffusion length (normalize first to 1)
hill = Inf;
noise = 0;

% initial conditions
p0 = [0.2 0.6];
iniON = round(p0*N);
I0 = [0 0];
dI = 0.01;
InitiateI = 0; % 0: no, 1: yes

% generate cell_type (0 case type 1, 1 case type 2)
cell_type = zeros(N,1);

% simulation parameters
tmax = 200;
mcsteps = 10^3;

% pos, dist
Lx = 1;
R = rcell*Lx/(gz+1); % disc radius
[pos, dist] = initial_cells_random_markov_periodic(...
    gz, Lx, R, mcsteps);

%{
fname_str = strrep(sprintf('N%d_iniON_%d_%d_M_int_%d_%d_%d_%d_a0_%.1f_Con_%d_%d_K_%d_%d_%d_%d_lambda_%.1f_%.1f', ...
    N, iniON(1), iniON(2), M_int(1,1), M_int(1,2), M_int(2,1), M_int(2,2), ...
    a0, Con(1), Con(2), K(1,1), K(1,2), K(2,1), K(2,2),...
    lambda(1), lambda(2)), '.', 'p');
%}

% check parameters
idx = (M_int == 0);
if ~all(K(idx)==0)
    fprintf('K has wrong entries! \n');
    warning('K has wrong entries!');
end

%}
%{
% TO DO: vectorize
dist_vec = a0*dist(1,:);
r = dist_vec(dist_vec>0); % exclude self influence
fN1 = sum(sinh(Rcell)*sum(exp((Rcell-r)./lambda(1)).*(lambda(1)./r)) ); % calculate signaling strength
fN2 = sum(sinh(Rcell)*sum(exp((Rcell-r)./lambda(2)).*(lambda(2)./r)) ); % calculate signaling strength

% nearest neighbour interaction strength
fprintf('activator fij(a0) = %.4f \n', sinh(Rcell)*sum(exp((Rcell-a0)./lambda(1)).*(lambda(1)./a0)))
fprintf('inhibitor fij(a0) = %.4f \n', sinh(Rcell)*sum(exp((Rcell-a0)./lambda(2)).*(lambda(2)./a0)))
%}

%% (2) Load parameters from saved trajectory
%
% with parameters saved as structure array 
% load data
%data_folder = 'H:\My Documents\Multicellular automaton\app\Multicellularity-2.1\data\time_evolution';
data_folder = 'D:\Multicellularity\app\git_repository\raw_current\data\time_evolution';
[file, path] = uigetfile(fullfile(data_folder, '\*.mat'), 'Load saved simulation');
load(fullfile(path, file));

s = save_consts_struct;
N = s.N;
a0 = s.a0;
K = s.K;
Con = s.Con;
Coff = s.Coff;
M_int = s.M_int;
hill = s.hill;
noise = s.noise;
rcell = s.rcell;
cells = cells_hist{1};
lambda = [1 s.lambda12];
%p0 = s.p_ini
tmax =  s.tmax;
gz = sqrt(N);
Rcell = rcell*a0;

% simulation parameters
%tmax = 10^4;
cell_type = zeros(N,1);
 
% Initial I
InitiateI = 0;
I0 = [0 0];
s_fields = fieldnames(s);
for i=1:numel(s_fields)
    if strcmp(s_fields{i},'I_ini_str')
        if ~isempty(s.I_ini_str)
            I0 = s.I_ini;
            InitiateI = 1;
        end
    end
end

% randomized pos, dist
%[dist, pos] = init_dist_hex(gz, gz);
Lx = 1;
R = rcell*Lx/(gz+1); % disc radius
mcsteps = 10^4;
[pos, dist] = initial_cells_random_markov_periodic(...
    gz, Lx, R, mcsteps);
%}
% Check whether loaded trajectory is same as simulation
%{
eq = 2*ones(numel(cells_hist_2), 1);
for i=1:numel(cells_hist_2)
    eq(i) = all(all(cells_hist{i} == cells_hist_2{i}));
end
%}
%% Simulate
% settings
disp_mol = 12;
showI = 0;

% generate initial lattice
%{
iniON = round(p0*N);
cells = zeros(N, 2);
for i=1:numel(iniON)
    cells(randperm(N,iniON(i)), i) = 1;
    if InitiateI && hill==Inf
        fprintf('Generating lattice with I%d(t=0)... \n', i);
        [cells_temp, test, I_ini] = generate_I_new(cells(:, i), I0(i), I0(i)+dI, dist, a0);
        cells(:,i) = cells_temp;
        fprintf('Generated initial I%d: %.2f; in range (Y=1/N=0)? %d; \n', i, I_ini, test);
    end
end
%}

% store initial config
cells_hist = {};
cells_hist{end+1} = cells; %{cells(:, 1), cells(:, 2)};

% dynamics
hin = figure();
t = 0;
plot_handle = reset_cell_figure(hin, pos, rcell);
update_figure_periodic_scatter(plot_handle, cells, t, disp_mol, showI, a0, dist);
[cellsOut, changed] = update_cells_two_signals_multiply_finite_Hill(cells, dist, M_int, a0,...
        Rcell, Con, Coff, K, lambda, hill, noise);
while changed && t < tmax
    pause(0.2);
    t = t+1;
    cells = cellsOut;
    cells_hist{end+1} = cells; %{cells(:, 1), cells(:, 2)};
    update_figure_periodic_scatter(plot_handle, cells, t, disp_mol, showI,...
        a0, dist);
    [cellsOut, changed] = update_cells_two_signals_multiply_finite_Hill(...
        cells, dist, M_int, a0, Rcell, Con, Coff, K, lambda, hill, noise);
end
t_out = t; % save final time

%% Plot p(t)
t0 = 0;
fig_pos = [1 1 10 8];
plot_p_vs_t(cells_hist, t0, fig_pos)

%% Plot I(t)
option = 1;
plot_I_vs_t(cells_hist, t0, a0, dist, option, fig_pos)