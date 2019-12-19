clear all
close all
%% Simulation parameters
% other settings
network = 15;
networks_all = [15 19 33 34 36];
tmax = 10^4; % max. number of time steps 

% default file name
sim_ID = 'two_signal_mult';

% (2) Load parameters that spontaneously generate TWs from batch simulations 
folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\batch_sim_all_topologies_run2';
fname_str = 'batch_sim_all_topologies_run2_count_TW_analyzed';

% load K, Con parameters
load(fullfile(folder, fname_str), 'Con_all_by_network', 'K_all_by_network');
network_idx = find(networks_all==network, 1);

Con_wave_sim = Con_all_by_network{network_idx};
K_wave_sim = K_all_by_network{network_idx};

Con_wave_sim = unique(Con_wave_sim, 'rows');
num_psets = size(Con_wave_sim, 1);
K_wave_sim = unique(K_wave_sim(:,:), 'rows');
K_wave_sim = reshape(K_wave_sim, num_psets, 2, 2);

fprintf('Number of parameter sets: %d \n',...
    num_psets);

% load other parameters
fname_str = 'batch_sim_analyzed_data_batch2';
load(fullfile(folder, fname_str), 'N', 'a0', 'hill', 'noise', 'lambda',...
    'InitiateI', 'rcell', 'M_int_all_reduced');

% choose random Con, K values
% Con = Con_wave_sim(1,:);
% K = squeeze(K_wave_sim(1,:,:));

lambda12 = lambda(2);
Coff = [1 1];
Rcell = rcell*a0;
gz = sqrt(N);
M_int = M_int_all_reduced{network};

%InitiateI = 0;
cell_type = zeros(N,1);

mcsteps = 0;
growth_rate = 0;
R_division = 0;
sigma_D = 0;

p0 = Inf;
I0 = Inf;
%% get all files in folder
for idx_param_loop = 2:10
    load_folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\trav_wave_with_noise\TW_formation_network_15';
    subfolder = sprintf('Parameter_set_%d', idx_param_loop);
    folder_orig_sim = fullfile(load_folder, subfolder);

    listing = dir(folder_orig_sim);
    num_files = numel(listing)-2;
    names_orig_sim = {};
    filecount = 0;
    %pattern = strrep(sprintf('%s_N%d_ini_state_TW_params_%d_t_out_%s_period_%s',...
    %        sim_ID, N, idx_param_loop, '(\d+)', '(\d+|Inf)' ), '.', 'p');
    %pattern = strrep(sprintf('%s_N%d_ini_state_TW_params_%d_sim_%s_t_out_%s_period_%s',...
    % 	sim_ID, N, idx_param_loop, '(\d+)', '(\d+)', '(\d+|Inf)' ), '.', 'p');
    
    pattern = strrep(sprintf('%s_N%d_params_%d_sim_%s_noise_%.3f_sim_%s_t_out_%s_period_%s',...
     	sim_ID, N, idx_param_loop, '(\d+)', 0, '(\d+)', '(\d+)', '(\d+|Inf)' ), '.', 'p');
    
    for i = 1:num_files
        filename = listing(i+2).name;
        % remove extension and do not include txt files
        [~,name,ext] = fileparts(filename);
        %if strcmp(ext, '.mat')
            match = regexp(name, pattern, 'match');
            %disp(match);
            if ~isempty(match)
                filecount = filecount + 1;
                names_orig_sim{end+1} = name;
            end
        %end
    end
    %%
    for i=1:numel(names_orig_sim)
        oldname = fullfile(folder_orig_sim, names_orig_sim{i});

        pattern = sprintf('%s_N%d_params_%d_sim_%s_noise_%.3f_sim_%s_%s',...
            sim_ID, N, idx_param_loop, '\d+', 0, '\d+', '(.+)' );

        [token, ~] = regexp(names_orig_sim{i}, pattern, 'tokens', 'match');
        %newname_str = sprintf('%s_N%d_ini_state_rand_params_%d_sim_%d_%s', sim_ID,...
        %    N, idx_param_loop, i, token{1}{1});
        newname_str = sprintf('%s_N%d_params_%d_sim_%d_noise_0p000_%s', sim_ID,...
            N, idx_param_loop, i, token{1}{1});

        newname = fullfile(folder_orig_sim, newname_str);
        disp(newname)
        movefile(strcat(oldname, '.mat'), strcat(newname, '.mat'));
    end

end
