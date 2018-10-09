%% Selects files of the simulation set of all topologies
clear variables
close all
clc

%% Parameters and settings
% Settings
% Note: increasing nsim at n_pset is always possible. However, increasing
% n_pset leads to data sets that do not form a perfect LHS sample
n_pset = 10000; % number of parameter sets to do
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

% get pos, dist
mcsteps = 0;
[pos, dist] = initial_cells_random_markov_periodic(gz, mcsteps, rcell);

% folders
%parent_folder = 'L:\BN\HY\Shared\Yiteng\two_signals\batch_sim_all_topologies';
parent_folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\batch_sim_all_topologies_run2';

%% Load trajectories, loop over all topologies
%{
M_int_all = {};
M_int_all_reduced = {};

%{
K_all = zeros(3^4, n_pset, 2, 2); % costly, could also make first dim. 45
Con_all = zeros(3^4, n_pset, 2); 

t_out_all = zeros(3^4, n_pset, nsim);
period_all = zeros(3^4, n_pset, nsim);
non_uniform_all = zeros(3^4, n_pset, nsim);
%}

% Troubleshooting
corrupt_files = {};
missing_file_count = zeros(3^4, n_pset);

done = zeros(3,3,3,3); % keeps track of which topologies have been found already (up to symmetry)
network_all = [];
for network=1:3^4 % network 1 = [0 0; 0 0]
    subfolder1 = sprintf('Network_%d', network);

    [i11, i12, i21, i22] = ind2sub([3, 3, 3, 3], network);
    M = [0 1 -1];
    M_int = [M(i11) M(i12); M(i21) M(i22)];
    M_int_all{end+1} = M_int;

    gM = [i22 i21; i12 i11];
    if done(i11,i12,i21,i22)
    	continue
    elseif network==1
        continue
    else
        M_int_all_reduced{end+1} = M_int;
        network_all(end+1) = network;
        done(i11,i12,i21,i22) = 1;
        done(gM(1,1),gM(1,2),gM(2,1),gM(2,2))=1;
    end
    
    %disp(network);
    %
    for idx1=1:n_pset
        subfolder2=sprintf('Param_%d', idx1);

        % get parameters
        fname = fullfile(parent_folder, subfolder1, subfolder2, 'parameters.mat');
        load(fname);
        disp(fname);

        K_all(network, idx1,:,:) = thisK;
        Con_all(network, idx1, :) = thisCon;

        for idx2=1:nsim
            % load file
            %disp(idx2);
            fname_str = sprintf('all_topologies_simulate-v%d.mat', idx2);
            fname = fullfile(parent_folder, subfolder1, subfolder2, fname_str);
            if exist(fname, 'file')~=2
                fname_str = sprintf('all_topologies_simulate-v%d_tmax_reached.mat', idx2);
                fname = fullfile(parent_folder, subfolder1, subfolder2, fname_str);
            end
            if exist(fname, 'file')~=2
                warning('Could not find file %s', fname);
                missing_file_count(network, idx1) = missing_file_count(k, idx1) + 1;
                continue
            end
            
            % try loading file
            try 
                load(fname);
            catch ME
                warning('Unable to load file %s', fname);
                %if strcmp(ME.identifier, 'MATLAB:load:unableToReadMatFile')
                corrupt_files{end+1} = fname; % store corrupt files for deletion
                continue
                %end
            end
            %disp(fname);
            
            % store variables
            t_out_all(network, idx1, idx2) = t_out;
            period_all(network, idx1, idx2) = period;
            non_uniform_all(network, idx1, idx2) = size(unique(cells_hist{end}, 'rows'), 1);
        end
    end
    %
end
n_networks = numel(M_int_all_reduced);
%}

%% Load saved data
folder = 'H:\My Documents\Multicellular automaton\data\two_signals\batch_sim_all_topologies';
fname_str = sprintf('batch_sim_analyzed_data_batch2.mat');
load(fullfile(folder, fname_str));
n_networks = numel(M_int_all_reduced);

%% Select data
%selection_idx = zeros(size(t_out_all));

% selection criteria
idx1 = t_out_all==tmax;
idx2 = period_all>4 & period_all<Inf;
selection_idx = idx1 | idx2;

%% Move files
parent_folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\batch_sim_all_topologies_run2\selected';
for i=1:numel(network_all)
    network = network_all(i);
    
    % create/select right subfolder
    subfolder = sprintf('Network_%d', network);
    move_folder = fullfile(parent_folder, subfolder);
    if exist(move_folder, 'dir')~=7
        fprintf('Create folder: %s \n', move_folder);
        mkdir(move_folder);
    end
    
    % move files
    
    
end

%% 
copyfile viridis.m move_folder