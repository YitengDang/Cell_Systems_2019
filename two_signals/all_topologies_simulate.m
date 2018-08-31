%% Obtain statistics on all possible topologies by simulation
clear variables
close all
clc

%% Loop over topologies and phases
% Settings
% single_cell = 0;
% sym = 0; % include symmetries? 0: symmetric diagrams are excluded, 1: everything included
nvar = 3; % number of variables to include
nsim = 10; %100; % number of simulations per parameter set

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

% 
parent_folder = 'D:\data\all_topologies_simulate';

% main output
%phases_all = {}; % store all phases, stored as matrix with phase of each interaction (1-6), 0 if no interaction (M_int(i,j)==0)
%state_diagrams = {}; % cell(Ns2, 1); % store all state diagrams (graph transition matrices)
%steady_states = {}; %cell(Ns2, 1); % store all steady states
%cycles_all = {}; %cell(Ns2, 1); % store all loop structures

%% test simulation
%{
M_int = [0 -1; 1 0];
K = [0 25; 15 0];
Con = [19 25];
Coff = [1 1];
parent_folder = 'D:\data\all_topologies_simulate';
p0 = [0.5 0.5];
InitiateI = 0;
[cells_hist, period, t_onset] = time_evolution_save_func(N, a0, Rcell, lambda,...
    hill, noise, M_int, K, Con, Coff, dist, p0, InitiateI, parent_folder);

% Test sim counter
max_trials = 100;
folder = 'D:\data\all_topologies_simulate';
pattern = 'all_topologies_simulate-v(\d+)';
[sim_todo, filecount] = all_topologies_simulate_count_todo(...
    max_trials, folder, pattern);
%}

%% loop over all phases
M = [0 1 -1]; % index to interaction
for k=38 %1:3^4
    disp(k);
    [i11, i12, i21, i22] = ind2sub([3, 3, 3, 3], k);
    % matrix associated with indices
    M_int = [M(i11) M(i12); M(i21) M(i22)];
    
    % ---- Generate random parameter set ----
    % bounds
    K_b = [1 40]; 
    Con_b = [1 40];
    
    % Latin hypercube
    idxK = find(M_int~=0);
    idxCon = find(sum(abs(M_int), 1)~=0);
    
    nK = numel(idxK); %sum(sum(abs(M_int))); 
    nCon = numel(idxCon); %sum(sum(abs(M_int), 1)>0);
    
    if nK+nCon==0
        continue
    end
    
    x = lhsdesign(nvar, nK+nCon);
    
    % Visualize parameters
    %{
    figure;
    hold on
    xlim(K_b);
    ylim(Con_b);
    %}
    for idx1=1:nvar
        thisK = zeros(2);
        thisCon = zeros(1,2);
        thisK(idxK) = (K_b(2) - K_b(1))*x(idx1, 1:nK) + K_b(1);
        thisCon(idxCon) = (Con_b(2) - Con_b(1))*x(idx1, nK+1:end) + Con_b(1); 
        
        % Visualize parameters
        % plot(thisK(2,2), thisCon(2), 'bo');
        % plot(thisK(1,2), thisK(2,1), 'bo');
        % plot(thisCon(1), thisCon(2), 'ro');
        
        % get save folder
        subfolder1 = sprintf('Network_%d', k);
        subfolder2 = sprintf('Param_%d', idx1);
        save_folder = fullfile(parent_folder, subfolder1, subfolder2);
        
        if exist(save_folder, 'dir')~=7
            mkdir(save_folder);
        end
        
        % save / load parameter set
        fname = fullfile(save_folder, 'parameters.mat');
        if exist(fname, 'file')==2
            load(fname);
        else
            save(fname, 'thisK', 'thisCon');
        end    
        
        % count # simulations to do
        pattern = 'all_topologies_simulate-v(\d+)';
        [sim_todo, ~] = all_topologies_simulate_count_todo(...
            nsim, save_folder, pattern);
        
        % simulate trajectory
        for count=1:sim_todo
            %disp(count);
            p0 = [0.5 0.5];
            [cells_hist, period, t_onset] = time_evolution_save_func(N, a0, Rcell, lambda,...
                hill, noise, M_int, thisK, thisCon, Coff, dist, p0, InitiateI, save_folder);
        end
    
    end
    
end
