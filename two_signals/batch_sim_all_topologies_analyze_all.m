%% Obtain statistics on all possible topologies by simulation
% Count fraction of travelling wave trajectories
% Based on pre-selection of trajectories
clear variables
close all
clc
set(0, 'defaulttextinterpreter', 'tex');
%% Parameters and settings
% Settings
remote = 0;

% Note: increasing nsim at n_pset is always possible. However, increasing
% n_pset leads to data sets that do not form a perfect LHS sample
n_pset = 10^4; % number of parameter sets to do
nsim = 10; % number of simulations per parameter set
tmax = 10000;

% Fixed parameters
gz = 15;
N = gz^2;
a0 = 1.5;
rcell = 0.2;
Rcell = rcell*a0;
lambda = [1 1.2];
hill = Inf;
noise = 0;
Coff = [1 1];
InitiateI = 0;
digits = -floor(log10(1/N)); % number of significant digits required to determine constancy of p

% get pos, dist
mcsteps = 0;
[pos, dist] = initial_cells_random_markov_periodic(gz, mcsteps, rcell);

% folders
parent_load_folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals_OR_logic\batch_sim_all_topologies';
save_folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals_OR_logic\batch_sim_all_topologies';
if remote
    parent_load_folder = strrep(parent_load_folder, 'N:\', 'W:\staff-bulk\');
    save_folder = strrep(save_folder, 'H:\', 'W:\staff-homes\d\yitengdang\');
end

% Get all networks
done = zeros(3,3,3,3); % keeps track of which topologies have been found already (up to symmetry)
network_all = [];
for network=1:3^4 % network 1 = [0 0; 0 0]
    subfolder1 = sprintf('Network_%d', network);
    [i11, i12, i21, i22] = ind2sub([3, 3, 3, 3], network);
    gM = [i22 i21; i12 i11];
    if done(i11,i12,i21,i22)
    	continue
    elseif network==1
        continue
    else
        network_all(end+1) = network;
        done(i11,i12,i21,i22) = 1;
        done(gM(1,1),gM(1,2),gM(2,1),gM(2,2))=1;
    end
end
%% Load full data
%networks_sel = [15 19 33 34 36]; %[15  16	19	20	32	33	34	36	43];
n_networks = 44;
t_out_all = zeros( numel(network_all), n_pset, nsim );
period_all = zeros( numel(network_all), n_pset, nsim );
non_uniform_all = zeros( numel(network_all), n_pset, nsim );
p_final_all = zeros( numel(network_all), n_pset, nsim, 2 );
I_final_all = zeros( numel(network_all), n_pset, nsim, 2 );
TW_count_strict = zeros( numel(network_all), n_pset );
TW_count_loose = zeros( numel(network_all), n_pset );
TW_K_by_network = cell(numel(network_all), 1);
TW_Con_by_network = cell(numel(network_all), 1);

done = zeros(3,3,3,3); % keeps track of which topologies have been found already (up to symmetry)
%network_all = [];
%%
for network_idx=15:numel(network_all) % network 1 = [0 0; 0 0]
    network = network_all(network_idx);
    subfolder1 = sprintf('Network_%d', network);
    %disp(network);
    
    K_found_TW = [];
    Con_found_TW = [];
    for idx1=1:n_pset
        subfolder2=sprintf('Param_%d', idx1);

        % get parameters
        fname = fullfile(parent_load_folder, subfolder1, subfolder2, 'parameters.mat');
        load(fname);
        disp(fname);

        K_all(network_idx, idx1, :, :) = thisK;
        Con_all(network_idx, idx1,  :) = thisCon;
        
        for idx2=1:nsim
            % load file
            %disp(idx2);
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
            
            % store variables
            %
            %disp(t_out);
            %disp(period);
            %pause(0.5);
            
            t_out_all(network, idx1, idx2) = t_out;
            period_all(network, idx1, idx2) = period;
            non_uniform_all(network, idx1, idx2) = size(unique(cells_hist{end}, 'rows'), 1);
            %}
            if period<Inf
                %disp('found');
                % average over one period
                p_avg = zeros(1,2);
                I_avg = zeros(1,2);
                for t1=numel(cells_hist):-1:numel(cells_hist)-period+1
                    p_avg = p_avg + mean(cells_hist{t1}, 1);
                    I_avg(1) = I_avg(1) + moranI(cells_hist{t1}(:,1), a0*dist);
                    I_avg(2) = I_avg(2) + moranI(cells_hist{t1}(:,2), a0*dist);
                end
                p_final_all(network, idx1, idx2, :) = p_avg/period;
                I_final_all(network, idx1, idx2, :) = I_avg/period;
            else
                % take final value
                p_final_all(network, idx1, idx2, 1) = mean(cells_hist{end}(:,1));
                p_final_all(network, idx1, idx2, 2) = mean(cells_hist{end}(:,2));
                I_final_all(network, idx1, idx2, 1) = moranI(cells_hist{end}(:,1), a0*dist);
                I_final_all(network, idx1, idx2, 2) = moranI(cells_hist{end}(:,2), a0*dist);
            end
            
            % Check if it is a TW, and store data if so
            if mod(period, gz)==0
                % Check whether selected files are TWs
                %disp('TW found!');
                %disp(filename);
                %a0 = save_consts_struct.a0;
                [trav_wave, trav_wave_2] = travelling_wave_test(cells_hist, a0,...
                    period, numel(cells_hist)-1, dist, digits);
                TW_count_strict(network_idx, idx2) = ...
                    TW_count_strict(network_idx, idx2) + trav_wave;
                TW_count_loose(network_idx, idx2) = ...
                    TW_count_loose(network_idx, idx2) + trav_wave_2;

                % Store K, Con values
                if trav_wave_2
                    count_TW = sum(TW_count_loose(network_idx, :));
                    K_found_TW(count_TW, :, :) = save_consts_struct.K;
                    Con_found_TW(count_TW, :) = save_consts_struct.Con;
                end
            end
            
        end % nsim loop
    end % pset loop
    % store found Con, K values
    TW_K_by_network{network_idx} = K_found_TW;
    TW_Con_by_network{network_idx} = Con_found_TW;
    
    %% Save data per network
    TW_count_loose_network = TW_count_loose(network_idx, :);
    TW_count_strict_network = TW_count_strict(network_idx, :);
    t_out_network = squeeze(t_out_all(network_idx, :, :));
    period_network = squeeze(period_all(network_idx, :, :));
    non_uniform_network = squeeze(non_uniform_all(network_idx, :, :));
    p_final_network = squeeze(p_final_all(network_idx, :, :, :));
    I_final_network = squeeze(I_final_all(network_idx, :, :, :));
    save_folder_data = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals_OR_logic\batch_sim_all_topologies\analyzed_data';
    fname_str = sprintf('batch_sim_OR_logic_all_topologies_analyzed_network_%d', network);
    save( fullfile(save_folder_data, fname_str), 'TW_count_strict_network',...
        'TW_count_loose_network', 'K_found_TW', 'Con_found_TW',...
        't_out_all', 'period_all', 'non_uniform_all', 'p_final_all', 'I_final_all');    
end % network loop
%}
%% Save all analyzed data
%
save_folder_data = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals_OR_logic\batch_sim_all_topologies';
fname_str = 'batch_sim_OR_logic_all_topologies_analyzed';
save( fullfile(save_folder_data, fname_str), 'TW_count_strict',...
    'TW_count_loose', 'K_all_by_network', 'Con_all_by_network',...
    't_out_all', 'period_all', 'non_uniform_all', 'p_final_all', 'I_final_all',...
    'TW_K_by_network', 'TW_Con_by_network' );
%}
%% Load analyzed data
%{
save_folder_data = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals_OR_logic\batch_sim_all_topologies';
fname_str = 'batch_sim_OR_logic_all_topologies_analyzed_by';
load( fullfile(save_folder_data, fname_str), 'networks_sel', 'TW_count_strict',...
    'TW_count_loose', 'K_all_by_network', 'Con_all_by_network',...
    't_out_all', 'period_all', 'non_uniform_all', 'p_final_all', 'I_final_all',...
    'TW_K_by_network', 'TW_Con_by_network' );
%}