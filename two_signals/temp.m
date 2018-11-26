clear all
close all
%%
N = 225;
K12_all = [5 8 9 15];

% Load file names
parent_folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\sweep K12 new lattice';
names = {};
subfolder = sprintf('N%d', N);
listing = dir(fullfile(parent_folder, subfolder));
num_files = numel(listing)-2; %first two entries are not useful
count = 0;
for i = 1:num_files
    filename = listing(i+2).name;
    % remove extension and do not include txt files
    [~,fname,ext] = fileparts(filename);
    if strcmp(ext, '.mat')
        count = count + 1;
        names{count} = fname;
    end
end

% File data on parameters
pattern = sprintf('two_signal_mult_N%d_initiateI0_randpos_mcsteps0_K12_%d_t_out_%s_period_%s-v%s',...
    N, K12, '(\d+)', '(\d+|Inf)', '(\d+)');
fnames_match = {};
for i=1:numel(names)
    [~, tokens] = regexp(names{i}, pattern, 'match', 'tokens');
    if ~isempty(tokens)
        disp(i);
        t_out = tokens{1}{1};
        period = tokens{1}{2};
        v = tokens{1}{3};
        
        fnames_match{end+1} = names{i};
        
    end
end

% Pick random files
n_files = 10; % number of files to pick
idx_pick = randperm( numel(fnames_match), n_files ); % indices of files to pick
fnames_new = {};
for i=1:n_files
    idx_file = idx_pick(i);
    [~, tokens] = regexp(fnames_match{idx_file}, pattern, 'match', 'tokens');
    t_out = tokens{1}{1};
    period = tokens{1}{2};
    v = tokens{1}{3};

    load(fullfile(parent_folder, subfolder, fnames_match{idx_file}));
    
    % Save as new file
    save_fname_str = sprintf('two_signal_mult_N%d_K12_%d_t_out_%d_period_%d-v%d',...
        N, K12, t_out, period, v);
    save_fname = fullfile(parent_folder, 'sensitivity_init_cond', save_fname_str);
    fnames_new{end+1} = save_fname;
    
    save(save_fname, 'cells_hist', 'distances', 'positions',...
        'save_consts_struct', 't_onset', 't_out');
end

%% Run new simulations with initial conditions slightly deviating from the earlier ones
% set base parameters
load(fnames_new{1}); % load first simulation
a0 = save_consts_struct.a0;
K = save_consts_struct.K;
Con = save_consts_struct.Con;
Coff = save_consts_struct.Coff;
M_int = save_consts_struct.M_int;
hill = save_consts_struct.hill;
noise = save_consts_struct.noise;
p_ini = save_consts_struct.p_ini;
rcell = save_consts_struct.rcell;
lambda12 = save_consts_struct.lambda12;
sim_ID = 'two_signal_mult';
mcsteps = 0;
InitiateI = 0;
Rcell = a0*rcell;
lambda = [1 lambda12];
tmax = 10^5; 

num_cells_changed = 1; % number of cells to change
n_trials = 2;
for i=1:n_files 
    fname = fnames_new{i};
    load(fname);
    cells_ini = cells_hist{1};

    for j=1:n_trials % number of trials for the same initial state
        % load base file
        
        % ----------- simulation ----------------------------------------------
        parent_folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\sweep K12 new lattice\';
        subfolder = 'sensitivity_init_cond';
        save_folder = fullfile(parent_folder, subfolder);
        fname_save = strcat(fname, sprintf('_%dcells_changed', num_cells_changed));

        [cells_hist, period, t_onset] = time_evolution_save_func_efficient_checks_change_cells(...
            N, a0, Rcell, lambda, hill, noise, M_int, K, Con, Coff,...
            distances, positions, sim_ID, mcsteps, tmax,...
            save_folder, fname_save, cells_ini, num_cells_changed);
        %----------------------------------------------------------------------
    end
end
