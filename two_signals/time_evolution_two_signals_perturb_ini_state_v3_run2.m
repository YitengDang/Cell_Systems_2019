%% Runs simulations by changing single cell states and simulating the system with one or two states changed
% Input a folder with simulations to perturb
clear all
close all
%%
% Parameters
N = 225;
K12_all = [8];
cells_to_change = 1;
%% Loop over parameters
for loop_idx = 1:numel(K12_all)
    K12 = K12_all(loop_idx);
    
    % save folder
    parent_folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\sweep K12 new lattice';
    subfolder = 'sensitivity_init_cond';
    subfolder2 = sprintf('K12_%d_%dcell_changed', K12, cells_to_change);
    save_folder = fullfile(parent_folder, subfolder, subfolder2);
    
    %% Select new set of simulations
    %{
    %----------------------------------------------------------------------
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
            fprintf('%s \n', names{i});
            fnames_match{end+1} = names{i};
        end
    end
    
    %
    % Pick random files from folder, and resave them    
    n_files = 10; % number of files to pick
    idx_pick = randperm( numel(fnames_match), n_files ); % indices of files to pick
    fnames_new = {};
    for i=1:n_files
        idx_file = idx_pick(i);
        [~, tokens] = regexp(fnames_match{idx_file}, pattern, 'match', 'tokens');
        t_out = tokens{1}{1};
        period = tokens{1}{2};
        version = str2double(tokens{1}{3});

        load(fullfile(parent_folder, subfolder, fnames_match{idx_file}));

        % Save as new file
        save_fname_str = sprintf('two_signal_mult_N%d_K12_%d_t_out_%d_period_%d-v%d',...
            N, K12, t_out, period, version);
        subfolder2 = sprintf('K12_%d', K12);
        save_fname = fullfile(parent_folder, 'sensitivity_init_cond',...
            subfolder2, 'originals', save_fname_str);
        fnames_new{end+1} = save_fname;
        
        fprintf('%s \n', save_fname);
        %{
        save(save_fname, 'cells_hist', 'distances', 'positions',...
            'save_consts_struct', 't_onset', 't_out');
        %}
    end
    %----------------------------------------------------------------------
    %}
    %% Load old simulations
    folder = fullfile(save_folder, 'originals');
    pattern = sprintf('two_signal_mult_N%d_K12_%d_t_out_%s_period_%s-v%s$',...
        N, K12, '(\d+)', '(\d+|Inf)', '(\d+)');
    
    fnames_new = {};
    listing = dir(folder);
    num_files = numel(listing)-2; %first two entries are not useful
    count = 0;
    for i = 1:num_files
        filename = listing(i+2).name;
        % remove extension and do not include txt files
        [~,fname_str,ext] = fileparts(filename);
        if strcmp(ext, '.mat')
            %[~, tokens] = regexp(fname pattern, 'match');
            if ~isempty(regexp(fname_str, pattern, 'match'))
                count = count + 1;
                fnames_new{count} = fullfile(folder, strcat(fname_str, '.mat'));
            end
        end
    end
    %% Count how many simulations to do
    n_trials = 200;
    sim_done = zeros(numel(fnames_new), 1);
    
    num_orig = numel(fnames_new); % number of original files
    
    % get base names
    fname_base_all = cell(numel(fnames_new), 1);
    for ii=1:num_orig
        %fname_base = fnames_new{ii};
        [~, fname_base_all{ii}] = fileparts(fnames_new{ii});
    end
    %
    listing = dir(save_folder);
    for i=1:numel(listing)-2
        fname = listing(i+2).name;
        for ii=1:num_orig
            sim_done(ii) = sim_done(ii) + ...
                ~isempty(regexp(fname,...
                fname_base_all{ii}, 'match'));
        end
    end
    %
    %{
    i = 2;
    ii = 2;
    fname = listing(i+2).name
    fname_base_all{ii}
    regexp(fname, fname_base_all{ii}, 'match')
    %}
    %
    % print result
    for ii=1:num_orig
    	fprintf('Original file %d, sims to do: %d \n', ii, n_trials-sim_done(ii) );
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
    if sum(strcmp( fields(save_consts_struct), 'p0' ))
        p_ini = save_consts_struct.p0;
    elseif sum(strcmp( fields(save_consts_struct), 'p_ini' ))
        p_ini = save_consts_struct.p_ini;
    end
    rcell = save_consts_struct.rcell;
    lambda12 = save_consts_struct.lambda12;
    sim_ID = 'two_signal_mult';
    mcsteps = 0;
    InitiateI = 0;
    Rcell = a0*rcell;
    lambda = [1 lambda12];
    tmax = 10^5;

    num_cells_changed = 1; % number of cells to change
    %n_trials = 180;
    for i=1:numel(fnames_new) 
        
        fname = fnames_new{i};
        load(fname);
        cells_ini = cells_hist{1};

        for j=1:n_trials % number of trials for the same initial state'
            fprintf('K12 = %d, file = %d, trial = %d/%d \n', K12, i, j, n_trials-sim_done(i) );
            if sim_done(i)+j > n_trials
                continue
            end
            % ----------- simulation ----------------------------------------------
               
            [~, fname_str] = fileparts(fname);
            fname_save = fullfile(save_folder, sprintf('%s_%dcells_changed',...
                fname_str, num_cells_changed));

            [~, ~, ~] = time_evolution_save_func_efficient_checks_change_cells(...
                N, a0, Rcell, lambda, hill, noise, M_int, K, Con, Coff,...
                distances, positions, sim_ID, mcsteps, tmax,...
                fname_save, cells_ini, num_cells_changed);
            %----------------------------------------------------------------------
        end
    end
    
end