close all
clear all
% maxNumCompThreads(4);
% warning off
%% (1) input parameters
% lattice parameters
gz = 15;
N = gz^2;
rcell = 0.2;
a0 = 1.5;
Rcell = 0.2*a0;

% circuit parameters 
M_int = [0 1; -1 -1];
Con = [18 16];
Coff = [1 1];
K = [0 8; 11 4];% K(i,j): sensitivity of type i to type j molecules
lambda = [1 1]; % diffusion length (normalize first to 1)
lambda12 = lambda(2)/lambda(1);
hill = Inf;
noise = 0;

% initial conditions
%p0 = [0.5 0.5];
%iniON = round(p0*N);
iniON = [51 61; 61 52];
I0 = [0 0];
dI = 0.01;
InitiateI = 0; % 0: no, 1: yes

% generate cell_type (0 case type 1, 1 case type 2)
cell_type = zeros(N,1);

% simulation parameters
tmax = 100;
nruns = 10;
mcsteps = 0;

%
fname_lbl = strrep(sprintf('N%d_iniON_%d_%d_M_int_%d_%d_%d_%d_a0_%.1f_Con_%d_%d_K_%d_%d_%d_%d_lambda_%.1f_%.1f_mcsteps_%d', ...
    N, iniON(1), iniON(2), M_int(1,1), M_int(1,2), M_int(2,1), M_int(2,2), ...
    a0, Con(1), Con(2), K(1,1), K(1,2), K(2,1), K(2,2),...
    lambda(1), lambda(2), mcsteps), '.', 'p');
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
%{
% with parameters saved as structure array 
% load data
data_folder = 'H:\My Documents\Multicellular automaton\app\Multicellularity-2.1\data\time_evolution';
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
%p0 = s.p_ini;
tmax =  s.tmax;
gz = sqrt(N);
Rcell = rcell*a0;
[dist, pos] = init_dist_hex(gz, gz);

% simulation parameters
%tmax = 100;
tmax = 10^4;
nruns = 10;
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

%}
% Check whether loaded trajectory is same as simulation
%{
eq = 2*ones(numel(cells_hist_2), 1);
for i=1:numel(cells_hist_2)
    eq(i) = all(all(cells_hist{i} == cells_hist_2{i}));
end
%}

%% ----------- simulation ------------------------------------
cells_hist = {};

% settings
disp_mol = 12;
showI = 0;

% generate initial lattice
%[dist, pos] = init_dist_hex(gz, gz);
nodisplay = 1;
[pos, dist] = initial_cells_random_markov_periodic(gz, mcsteps, rcell, nodisplay);

% generate initial state
%iniON = round(p0*N);
cells = zeros(N, 2);
for i=1:numel(iniON)
    cells(randperm(N,iniON(i)), i) = 1;
    if InitiateI && hill==Inf
        %fprintf('Generating lattice with I%d(t=0)... \n', i);
        dI = 0.1;
        [cells_temp, test, I_ini] = generate_I_new(cells(:, i), I0(i), I0(i)+dI, dist, a0);
        cells(:,i) = cells_temp;
        %fprintf('Generated initial I%d: %.2f; in range (Y=1/N=0)? %d; \n', i, I_ini, test);
    end
end

% store initial config
cells_hist{end+1} = cells; %{cells(:, 1), cells(:, 2)};

%-------------dynamics-----------------------------------------
hin=figure;
plot_handle = reset_cell_figure(hin, pos, rcell);
t = 0;
period = Inf; %default values
t_onset = Inf; 
[cellsOut, changed] = update_cells_two_signals_multiply_finite_Hill(cells, dist, M_int, a0,...
        Rcell, Con, Coff, K, lambda, hill, noise);
update_figure_periodic_scatter(plot_handle, cells, t, disp_mol, showI, a0, dist);

% always check within first t_ac time steps
t_ac = 10^2; 
while changed && period==Inf && t<t_ac
    %disp(t);
    pause(1);
    t = t+1;
    cells = cellsOut;
    cells_hist{end+1} = cells; %{cells(:, 1), cells(:, 2)};
    [period, t_onset] = periodicity_test_short(cells_hist); 
    update_figure_periodic_scatter(plot_handle, cells, t, disp_mol, showI, a0, dist);
    [cellsOut, changed] = update_cells_two_signals_multiply_finite_Hill(cells, dist, M_int, a0,...
        Rcell, Con, Coff, K, lambda, hill, noise);
end
% check periodically after t_ac time steps, with period t_check
t_check = 10^3; 
while changed && period==Inf && t<tmax
    %disp(t);
    pause(0.2);
    t = t+1;
    cells = cellsOut;
    cells_hist{end+1} = cells; %{cells(:, 1), cells(:, 2)};
    if mod(t, t_check)==0
        [period, t_onset] = periodicity_test_short(cells_hist); 
    end
    %update_cell_figure_continuum(app, pos, dist, a0, cells, app.Time, cell_type, disp_mol, 0);
    [cellsOut, changed] = update_cells_two_signals_multiply_finite_Hill(cells, dist, M_int, a0,...
        Rcell, Con, Coff, K, lambda, hill, noise);
end
%pause(1);

t_out = t; %default t_out
% if periodicity found, refine check to find period
if period<Inf && t>t_ac
    [period, t_onset] = periodicity_test_detailed(cells_hist, t_check, period);
    t_out = t_onset + period; 
end

if changed && t==tmax
    tmax_string = '_tmax_reached';
else
    tmax_string = '';
end
fprintf('Final: t_out = %d, period %d \n', t_out, period);

%% Save trajectory
save_folder = 'D:\temp';

% default file name
I_ini_str = '';
if InitiateI
    I_ini_str = sprintf('_I_ini_%.2f_%.2f', I0(1), I0(2));
end
fname_str = strrep(sprintf('%s%s_t_out_%d_period_%s%s',...
    fname_lbl, I_ini_str, t_out, num2str(period), tmax_string), '.', 'p');
ext = '.mat';
label = '';
        
% check if filename already exists
i=1;
fname = fullfile(save_folder, strcat(fname_str, '-v', num2str(i), label, ext));
while exist(fname, 'file') == 2
    i=i+1;
    fname = fullfile(save_folder, strcat(fname_str, '-v', num2str(i), label, ext));
end

if InitiateI
    save_vars = {N, a0, K, Con, Coff, M_int, hill, noise, p0, I0, rcell,...
        lambda12, I_ini_str, mcsteps};
    save_vars_lbl = {'N', 'a0', 'K', 'Con', 'Coff', 'M_int', 'hill', 'noise', 'p_ini', 'I_ini', 'rcell',...
        'lambda12','I_ini_str', 'mcsteps'};
else
    save_vars = {N, a0, K, Con, Coff, M_int, hill, noise, p0, rcell,...
        lambda12, I_ini_str, mcsteps};
    save_vars_lbl = {'N', 'a0', 'K', 'Con', 'Coff', 'M_int', 'hill', 'noise', 'p_ini', 'rcell',...
        'lambda12', 'I_ini_str', 'mcsteps'};
end

save_consts_struct = cell2struct(save_vars, save_vars_lbl, 2);
positions = pos;
distances = dist;

qsave = 0;
if qsave
    save(fname, 'save_consts_struct', 'cells_hist', 't_out',...
        'changed', 'period', 't_onset', 'positions', 'distances');
        fprintf('Saved simulation: %s ; \n', fname);
end