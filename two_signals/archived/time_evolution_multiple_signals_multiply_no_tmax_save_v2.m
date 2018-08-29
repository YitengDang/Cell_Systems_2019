close all
clear all
maxNumCompThreads(3);
%warning off
%% Simulation parameters
max_trials = 500;

% lattice parameters
%gz = 15;
%N = gz^2;
gz_all = 15;

% loop over K12
K12_all = 4:28;
%% (2) Load parameters from saved trajectory
%{
% with parameters saved as structure array 
% load data
%data_folder = 'H:\My Documents\Multicellular automaton\app\git_repository\release_2_1_full\data\time_evolution\sample_trajectories\typeIV';
data_folder = 'H:\My Documents\Multicellular automaton\app\git_repository\raw_current\data\time_evolution\parameter set 2b';
%data_folder = 'D:\Multicellularity\app\Multicellularity-2.1\data\time_evolution';
%file = 'two_signal_mult_M_int1_1_-1_-1_chaotic_state_tmax5000-v1.mat';
%file = 'two_signal_mult_M_int0_1_-1_1_long_t_trans_wave_to_period10_trav_wave-v1';
%file = 'two_signal_mult_M_int0_1_-1_1_transient_wave_to_travelling_wave-v1';
file = 'two_signal_mult_N225_initiateI1_I_ini_0p50_0p50_t_out_2436_period_15-v1';
%[file, data_folder] = uigetfile(fullfile(data_folder, '\*.mat'), 'Load saved simulation');
load(fullfile(data_folder, file));

s = save_consts_struct;
%N = s.N;
a0 = s.a0;
K = s.K;
Con = s.Con;
Coff = s.Coff;
M_int = s.M_int;
hill = s.hill;
noise = s.noise;
rcell = s.rcell;
cells = cells_hist{1};
lambda12 = s.lambda12;
lambda = [1 lambda12];
%p0 = s.p_ini;
p0 = [0.5 0.5];
%tmax =  s.tmax;
%gz = sqrt(N);
Rcell = rcell*a0;

% simulation parameters
%tmax = 10^4;


% Initial I
InitiateI = 0;
%{
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
%}
%% Specify parameters by hand 

a0 = 1.5;
K = [0 4; 11 4];
Con = [18 16];
Coff = [1 1];
M_int = [0 1; -1 1];
hill = Inf;
noise = 0;
rcell = 0.2;
%cells = cells_hist{1};
lambda12 = 1.2;
lambda = [1 lambda12];
%p0 = s.p_ini;
p0 = [0.5 0.5];
%tmax =  s.tmax;
%gz = sqrt(N);
Rcell = rcell*a0;
InitiateI = 0;
%%
for gz_idx=1:numel(gz_all)    
    
gz = gz_all(gz_idx);
N = gz^2;

[dist, pos] = init_dist_hex(gz, gz);
cell_type = zeros(N,1);

for K12_idx=1:numel(K12_all)

K(1,2) = K12_all(K12_idx);
%% input parameters
%{
a0 = 1.5;
rcell = 0.2;
Rcell = rcell*a0;

% circuit parameters
Con = [18 16];
Coff = [1 1];
M_int = [1 1; -1 1];
K = [3 12; 13 20]; % K(i,j): sensitivity of type i to type j molecules
lambda = [1 1.2]; % diffusion length (normalize first to 1)
lambda12 = lambda(2)/lambda(1);
hill = Inf;
noise = 0;

% initial conditions
p0 = [0.6 0.4];
iniON = round(p0*N);
I0 = [0 0];
dI = 0.01;
InitiateI = 0; % 0: no, 1: yes

% generate cell_type (0 case type 1, 1 case type 2)
cell_type = zeros(N,1);

% simulation parameters
tmax = 10000;
%nruns = 80;

% pos, dist
[dist, pos] = init_dist_hex(gz, gz);
%[pos,ex,ey] = init_cellpos_hex(gridsize,gridsize);
%dist = dist_mat(pos,gridsize,gridsize,ex,ey);

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

I_ini_str = '';
if InitiateI
    I_ini_str = sprintf('_I1_I2_%.2f_%.2f', I0(1), I0(2));
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
%% Check existing files
% Count how many simulations have already been done
%folder = fullfile('L:\HY\Shared\Yiteng\two_signals\parameter set 2b', sprintf('N%d', N));
%folder = fullfile('L:\BN\HY\Shared\Yiteng\two_signals', 'sweep K12');
folder = fullfile('L:\HY\Shared\Yiteng\two_signals', 'sweep K12');

if exist(folder, 'dir') ~= 7
    warning('Folder does not exist! ');
end
sim_ID = 'two_signal_mult';

% default file name
I_ini_str = '';
if InitiateI
    I_ini_str = sprintf('_I_ini_%.2f_%.2f', I0(1), I0(2));
end

filecount = 0;
%pattern = strrep(sprintf('%s_N%d_initiateI%d_%s_t_out_%s_period_%s',...
%        sim_ID, N, InitiateI, I_ini_str, '(\d+)', '(\d+|Inf)'), '.', 'p');
pattern = strrep(sprintf('%s_N%d_initiateI%d%s_K12_%d_t_out_%s_period_%s',...
        sim_ID, N, InitiateI, I_ini_str, K(1,2), '(\d+)', '(\d+|Inf)'), '.', 'p');


listing = dir(folder);
num_files = numel(listing)-2;
names = {};
for i = 1:num_files
    filename = listing(i+2).name;
    % remove extension and do not include txt files
    [~,name,ext] = fileparts(filename);
    if strcmp(ext, '.mat')
        match = regexp(name, pattern, 'match');
        %disp(match);
        if ~isempty(match)
            filecount = filecount + 1;
            names{end+1} = name;
        end
    end
end

fprintf('N=%d, K12 = %d, sim to do: %d \n', N, K(1,2), max_trials-filecount);
%% Simulate

%{
folder = fullfile('L:\HY\Shared\Yiteng\two_signals');

sim_ID = 'two_signal_mult';
M_int_str = sprintf('M_int%d_%d_%d_%d', M_int(1,1), M_int(1,2), M_int(2,1), M_int(2,2));
K_str = sprintf('K%.1f_%.1f_%.1f_%.1f', K(1,1), K(1,2), K(2,1), K(2,2));
Con_str = sprintf('Con%.1f_%.1f', Con(1), Con(2));
Coff_str = sprintf('Coff%.1f_%.1f', Coff(1), Coff(2));

lambda12 = lambda(2);
p_ini_str = sprintf('p_ini%.2f_%.2f', p0(1), p0(2));
fname_str = strrep(sprintf('%s_%s_hill%.2f_N%d_a0_%.2f_%s_%s_%s_noise%.1f_%s%s_l12_%.1f_rcell%.1f_tmax%d',...
    sim_ID, M_int_str, hill, N, a0, K_str, Con_str, Coff_str, noise,...
    p_ini_str, I_ini_str, lambda12, rcell, tmax), '.', 'p');
%}

for trial=1:max_trials-filecount
    fprintf('trial %d \n', trial);
    cells_hist = {};
    %t_out = 0;
    %changed = 1;

    % generate initial lattice
    iniON = round(p0*N);
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
    t = 0;
    period = Inf; %default values
    t_onset = Inf; 
    [cellsOut, changed] = update_cells_two_signals_multiply_finite_Hill(cells, dist, M_int, a0,...
            Rcell, Con, Coff, K, lambda, hill, noise);

    while changed && period==Inf %&& t<tmax
        %pause(0.2);
        t = t+1;
        cells = cellsOut;
        cells_hist{end+1} = cells; %{cells(:, 1), cells(:, 2)};
        [period, t_onset] = periodicity_test_short(cells_hist); 

        %update_cell_figure_continuum(app, pos, dist, a0, cells, app.Time, cell_type, disp_mol, 0);
        [cellsOut, changed] = update_cells_two_signals_multiply_finite_Hill(cells, dist, M_int, a0,...
            Rcell, Con, Coff, K, lambda, hill, noise);
    end
    t_out = t; % save final time
    
    %if t_out < 4*2^N
    %    periodic = 'chaotic';
    %    fprintf('t_out = %d, no periodicity found \n', t_out);
    %else
    %    periodic = 'periodic';
    %    fprintf('t_out = %d, period %d \n', t_out, period);
    %end
    fprintf('t_out = %d, period %d \n', t_out, period);
    
    % temp
    %cells_hist = zeros(N,1);
    %t_onset = Inf;
    %period = Inf;
    %t_out = tmax;
    %periodic = 'periodic';
    %--------------------------------------------------------------
    % Save result
    fname_str = strrep(sprintf('%s_N%d_initiateI%d%s_K12_%d_t_out_%d_period_%s',...
        sim_ID, N, InitiateI, I_ini_str, K(1,2), t_out, num2str(period)), '.', 'p');
    ext = '.mat';
    label = '';

    % filename 
    %filename = strrep(sprintf('%s_%s_hill%.2f_N%d_a0_%.2f_%s_%s_%s_noise%.1f_%s%s_l12_%.1f_rcell%.1f_tmax%d',...
    %    sim_ID, M_int_str, hill, N, a0, K_str, Con_str, Coff_str, noise, p_ini_str, I_ini_str, lambda12, rcell, tmax), '.', 'p');

    % check if filename already exists
    i=1;
    fname = fullfile(folder, strcat(fname_str, '-v', num2str(i), label, ext));
    while exist(fname, 'file') == 2
        i=i+1;
        fname = fullfile(folder, strcat(fname_str, '-v', num2str(i), label, ext));
    end
    
    if InitiateI
        save_vars = {N, a0, K, Con, Coff, M_int, hill, noise, p0, I0, rcell,...
            lambda12, sim_ID, I_ini_str};
        save_vars_lbl = {'N', 'a0', 'K', 'Con', 'Coff', 'M_int', 'hill', 'noise', 'p_ini', 'I_ini', 'rcell',...
            'lambda12', 'sim_ID', 'I_ini_str'};
    else
        save_vars = {N, a0, K, Con, Coff, M_int, hill, noise, p0, rcell,...
            lambda12, sim_ID, I_ini_str};
        save_vars_lbl = {'N', 'a0', 'K', 'Con', 'Coff', 'M_int', 'hill', 'noise', 'p_ini', 'rcell',...
            'lambda12', 'sim_ID', 'I_ini_str'};
    end
    
    save_consts_struct = cell2struct(save_vars, save_vars_lbl, 2);

    save(fname, 'save_consts_struct', 'cells_hist', 't_out',...
        'changed', 'period', 't_onset');
        fprintf('Saved simulation: %s ; \n', fname);
    %--------------------------------------------------------------------------
end
%}
end
end