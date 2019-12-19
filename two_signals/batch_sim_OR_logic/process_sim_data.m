% Processes and saves the simulated data
clear all
close all
clc

%% Get all networks and M_int
done = zeros(3,3,3,3); % keeps track of which topologies have been found already (up to symmetry)
network_all = [];
M_int_all_reduced = {};
for network=1:3^4 % network 1 = [0 0; 0 0]
    subfolder1 = sprintf('Network_%d', network);
    [i11, i12, i21, i22] = ind2sub([3, 3, 3, 3], network);
    M = [0 1 -1];
    M_int = [M(i11) M(i12); M(i21) M(i22)];
    gM = [i22 i21; i12 i11];
    if done(i11,i12,i21,i22)
    	continue
    elseif network==1
        continue
    else
        network_all(end+1) = network;
        M_int_all_reduced{end+1} = M_int;
        done(i11,i12,i21,i22) = 1;
        done(gM(1,1),gM(1,2),gM(2,1),gM(2,2))=1;
    end
end

%% Setup
n_pset = 10000;
n_sim = 10;
n_network = numel(network_all);
parent_load_folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals_batch_sim_2\batch_sim_all_topologies_OR_logic_run1';

%% Load full data
K_all_network = zeros(n_network, n_pset, 2, 2);
Con_all_network = zeros(n_network, n_pset, 2);
I_final_all_network = zeros(n_network, n_pset, n_sim, 2);
p_final_all_network = zeros(n_network, n_pset, n_sim, 2);
non_uniform_all_network = zeros(n_network, n_pset, n_sim);
period_all_network = zeros(n_network, n_pset, n_sim);
t_out_all_network = zeros(n_network, n_pset, n_sim);
TW_Con_by_network = cell(n_network, 1);
TW_K_by_network = cell(n_network, 1);
TW_count_loose_network_all = zeros(n_network, n_pset);
TW_count_strict_network_all = zeros(n_network, n_pset);
%%
for network_idx=1:n_network
    network = network_all(network_idx);
    subfolder1 = sprintf('Network_%d', network);
    %disp(network);
    
    K_all = zeros( n_pset, 2, 2 );
    Con_all = zeros( n_pset, 2 );
    t_out_all = zeros( n_pset, n_sim );
    period_all = zeros( n_pset, n_sim );
    non_uniform_all = zeros( n_pset, n_sim );
    p_final_all = zeros( n_pset, n_sim, 2 );
    I_final_all = zeros( n_pset, n_sim, 2 );
    TW_count_strict = zeros( 1, n_pset );
    TW_count_loose = zeros( 1, n_pset );
    K_found_TW = [];
    Con_found_TW = [];
    
    fpattern = sprintf('batch_sim_OR_logic_network_%s_t_out_%s_period_%s-v%s',...
        '(\d+)', '(\d+)', '(\d+|Inf)', '(\d+)');
    for idx1=1:n_pset
        subfolder2=sprintf('Param_%d', idx1);
        
        % get parameters
        fname = fullfile(parent_load_folder, subfolder1, subfolder2, 'parameters.mat');
        load(fname);
        disp(fname);
        
        %K_all_network(network_idx, idx1, :, :) = thisK;
        %Con_all_network(network_idx, idx1,  :) = thisCon;
        K_all(idx1, :, :) = thisK;
        Con_all(idx1, :) = thisCon;
        
        % all all filenames
        foldercont = dir( fullfile(parent_load_folder, subfolder1, subfolder2) );
        names = {foldercont(3:end-1).name};
        
        for idx2=1:n_sim
            % load file
            load_fname = names{idx2};
            %disp(idx2);
            out = regexp(load_fname, fpattern, 'match');
            
            if ~isempty(out)
                fname = fullfile(parent_load_folder, subfolder1, subfolder2, load_fname);
                % try loading file
                try 
                    load(fname);
                    %disp(fname);
                catch ME
                    warning('Unable to load file %s', fname);
                    %if strcmp(ME.identifier, 'MATLAB:load:unableToReadMatFile')
                    corrupt_files{end+1} = fname; % store corrupt files for deletion
                    continue
                    %end
                end
                %disp(fname);
            else
                disp('Unable to find file');
            end
           
            %{
            fname_str = sprintf('all_topologies_simulate-v%d.mat', idx2);
            fname = fullfile(parent_load_folder, subfolder1, subfolder2, fname_str);
            if exist(fname, 'file')~=2
                fname_str = sprintf('all_topologies_simulate-v%d_tmax_reached.mat', idx2);
                fname = fullfile(parent_load_folder, subfolder1, subfolder2, fname_str);
            end
            if exist(fname, 'file')~=2
                warning('Could not find file %s', fname);
                missing_file_count(network, idx1) = missing_file_count(network, idx1) + 1;
                continue
            end
            %}
            a0 = save_consts_struct.a0;
            gz = sqrt(save_consts_struct.N);
            t_out_all(idx1, idx2) = t_out;
            period_all(idx1, idx2) = period;
            non_uniform_all(idx1, idx2) = size(unique(cells_hist{end}, 'rows'), 1);
            %}
            if period<Inf
                %disp('found');
                % average over one period
                p_avg = zeros(1,2);
                I_avg = zeros(1,2);
                for t1=numel(cells_hist):-1:numel(cells_hist)-period+1
                    p_avg = p_avg + mean(cells_hist{t1}, 1);
                    I_avg(1) = I_avg(1) + moranI(cells_hist{t1}(:,1), a0*distances);
                    I_avg(2) = I_avg(2) + moranI(cells_hist{t1}(:,2), a0*distances);
                end
                p_final_all(idx1, idx2, :) = p_avg/period;
                I_final_all(idx1, idx2, :) = I_avg/period;
            else
                % take final value
                p_final_all(idx1, idx2, 1) = mean(cells_hist{end}(:,1));
                p_final_all(idx1, idx2, 2) = mean(cells_hist{end}(:,2));
                I_final_all(idx1, idx2, 1) = moranI(cells_hist{end}(:,1), a0*distances);
                I_final_all(idx1, idx2, 2) = moranI(cells_hist{end}(:,2), a0*distances);
            end
            
            % Check if it is a TW, and store data if so
            if mod(period, gz)==0
                % Check whether selected files are TWs
                %disp('TW found!');
                %disp(filename);
                %a0 = save_consts_struct.a0;
                [trav_wave, trav_wave_2] = travelling_wave_test(cells_hist, a0,...
                    period, numel(cells_hist)-1, distances, digits);
                TW_count_strict(idx1) = ...
                    TW_count_strict(idx1) + trav_wave;
                TW_count_loose(idx1) = ...
                    TW_count_loose(idx1) + trav_wave_2;

                % Store K, Con values
                if trav_wave_2
                    count_TW = sum(TW_count_loose);
                    K_found_TW(count_TW, :, :) = save_consts_struct.K;
                    Con_found_TW(count_TW, :) = save_consts_struct.Con;
                end
            end
            
        end % nsim loop
    end % pset loop
    
    % store found Con, K values
    TW_K_by_network{network_idx} = K_found_TW;
    TW_Con_by_network{network_idx} = Con_found_TW;
    % store other variables
    K_all_network(network_idx, :, :, :) = K_all;
    Con_all_network(network_idx, :,  :) = Con_all;
    TW_count_loose_network_all(network_idx, :) = TW_count_loose;
    TW_count_strict_network_all(network_idx, :) = TW_count_strict;
    t_out_all_network(network_idx, :, :) = t_out_all(:, :);
    period_all_network(network_idx, :, :) = period_all(:, :);
    non_uniform_all_network(network_idx, :, :) = non_uniform_all(:, :);
    p_final_all_network(network_idx, :, :, :) = p_final_all(:, :, :);
    I_final_all_network(network_idx, :, :, :) = I_final_all(:, :, :);
    
    %% Save data per network
    save_folder_data = fullfile(parent_load_folder, 'analyzed_data');
    fname_str = sprintf('batch_sim_OR_logic_all_topologies_analyzed_network_%d', network);
    save( fullfile(save_folder_data, fname_str), 'K_all', 'Con_all',...
        'TW_count_strict', 'TW_count_loose', 'K_found_TW', 'Con_found_TW',...
        't_out_all', 'period_all', 'non_uniform_all', 'p_final_all', 'I_final_all');    
end % network loop
%}

%% Load and re-save all analyzed data

% Define variables
%{
K_all_network = zeros(n_network, n_pset, 2, 2);
Con_all_network = zeros(n_network, n_pset, 2);
I_final_all_network = zeros(n_network, n_pset, n_sim, 2);
p_final_all_network = zeros(n_network, n_pset, n_sim, 2);
non_uniform_all_network = zeros(n_network, n_pset, n_sim);
period_all_network = zeros(n_network, n_pset, n_sim);
t_out_all_network = zeros(n_network, n_pset, n_sim);
TW_Con_by_network = cell(n_network, 1);
TW_K_by_network = cell(n_network, 1);
TW_count_loose_network_all = zeros(n_network, n_pset);
TW_count_strict_network_all = zeros(n_network, n_pset);
%}
%{
for network_idx = 1 %:n_network
    network = network_all(network_idx);
    disp(network);
    
    % Load data
    load_folder_data = fullfile(parent_load_folder, 'analyzed_data');
    fname_str = sprintf('batch_sim_OR_logic_all_topologies_analyzed_network_%d', network);
    load( fullfile(load_folder_data, fname_str), ...
        'TW_count_strict', 'TW_count_loose', 'K_found_TW', 'Con_found_TW',...
        't_out_all', 'period_all', 'non_uniform_all', 'p_final_all', 'I_final_all');   
    
    % store found Con, K values
    TW_K_by_network{network_idx} = K_found_TW;
    TW_Con_by_network{network_idx} = Con_found_TW;
    % store other variables
    TW_count_loose_network_all(network_idx, :) = TW_count_loose;
    TW_count_strict_network_all(network_idx, :) = TW_count_strict;
    t_out_all_network(network_idx, :, :) = t_out_all(:, :);
    period_all_network(network_idx, :, :) = period_all(:, :);
    non_uniform_all_network(network_idx, :, :) = non_uniform_all(:, :);
    p_final_all_network(network_idx, :, :, :) = p_final_all(:, :, :);
    I_final_all_network(network_idx, :, :, :) = I_final_all(:, :, :);
end 
%}

% Resave data
%
save_folder_data = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals_batch_sim_2\batch_sim_all_topologies_OR_logic_run1\analyzed_data';
fname_str = 'batch_sim_OR_logic_all_topologies_analyzed_all';
save( fullfile(save_folder_data, fname_str), 'network_all', 'M_int_all_reduced',...
    'TW_count_strict_network_all', 'TW_count_loose_network_all',...
    'K_all_network', 'Con_all_network', 't_out_all_network', ...
    'period_all_network', 'non_uniform_all_network',...
    'p_final_all_network', 'I_final_all_network',...
    'TW_K_by_network', 'TW_Con_by_network' );
%}

%{
for network_idx = 1:numel(network_all)
    network = network_all(network_idx);
    fprintf('Network %d \n', network);
    
    %network = 2;
    load_folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals_batch_sim_2\batch_sim_all_topologies_OR_logic_run1'; 
    %'N:\tnw\BN\HY\Shared\Yiteng\two_signals_OR_logic\batch_sim_all_topologies\analyzed_data_old';
    %folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals_OR_logic\batch_sim_all_topologies\analyzed_data';
    fname_str = sprintf('batch_sim_OR_logic_all_topologies_analyzed_network_%d', network);
    load( fullfile(load_folder, fname_str) )

    %%
    I_final_all = squeeze(I_final_all(network, :, :, :));
    p_final_all = squeeze(p_final_all(network, :, :, :));
    non_uniform_all = squeeze(non_uniform_all(network, :, :));
    period_all = squeeze(period_all(network, :, :));
    t_out_all = squeeze(t_out_all(network, :, :));
    
    %{
    disp( p_final_all(1:20) );
    disp( I_final_all(1:20) );
    disp( t_out_all(1:20) );
    pause(1);
    
    save_folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals_OR_logic\batch_sim_all_topologies\analyzed_data';
    fname_out = fullfile(save_folder, fname_str);
    save( fname_out, 'I_final_all', 'p_final_all', 'non_uniform_all',...
        'period_all', 't_out_all', 'Con_found_TW', 'K_found_TW',...
        'TW_count_loose_network', 'TW_count_strict_network');
    %}
    %%
    I_final_all_network(network_idx, :, :, :) = I_final_all;
    p_final_all_network(network_idx, :, :, :) = p_final_all;
    non_uniform_all_network(network_idx, :, :, :) = non_uniform_all;
    period_all_network(network_idx, :, :, :) = period_all;
    t_out_all_network(network_idx, :, :, :) = t_out_all;
    Con_found_TW_network{network_idx} = Con_found_TW;
    K_found_TW_network{network_idx} = Con_found_TW;
    TW_count_loose_network_all(network_idx, :) = TW_count_loose_network;
    TW_count_strict_network_all(network_idx, :) = TW_count_strict_network;
end
%% 
save_folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals_OR_logic\batch_sim_all_topologies\analyzed_data';
fname_str = 'batch_sim_OR_logic_all_topologies_analyzed_all';
fname_out = fullfile(save_folder, fname_str);
%load(fname_out);

%%
save( fname_out, 'network_all', 'M_int_all_reduced', 'I_final_all_network', 'p_final_all_network', 'non_uniform_all_network',...
        'period_all_network', 't_out_all_network', 'Con_found_TW_network', 'K_found_TW_network');
%} 