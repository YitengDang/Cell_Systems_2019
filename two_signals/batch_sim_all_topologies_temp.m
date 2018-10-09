%% Obtain statistics on all possible topologies by simulation
% Loop over topologies and variables
clear variables
close all
clc

%% Settings
% Note: increasing nsim at n_pset is always possible. However, increasing
% n_pset leads to data sets that do not form a perfect LHS sample
n_pset = 10000; % number of parameter sets to do
nsim = 0; % number of simulations per parameter set

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

% Initial conditions
p0 = Inf; % special setting: random initial conditions
%p0 = [0.5 0.5];
I0 = Inf; 
tmax = 10000;
InitiateI = 0;

% get pos, dist
mcsteps = 0;
[pos, dist] = initial_cells_random_markov_periodic(gz, mcsteps, rcell);

%parent_folder = 'L:\BN\HY\Shared\Yiteng\two_signals\batch_sim_all_topologies';
parent_folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\batch_sim_all_topologies_run2';

%% test simulation
%{
M_int = [0 -1; 1 0];
K = [0 25; 15 0];
Con = [19 25];
Coff = [1 1];
%parent_folder = 'D:\data\all_topologies_simulate';
parent_folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\batch_sim_all_topologies_run2';
InitiateI = 0;
tmax = 10000;
[cells_hist, period, t_onset] = time_evolution_save_func(...
    N, a0, Rcell, lambda, hill, noise, M_int, K, Con, Coff,...
    dist, InitiateI, p0, I0, tmax, parent_folder);
%% Test sim counter
max_trials = 100;
folder = 'D:\data\all_topologies_simulate';
pattern = 'all_topologies_simulate-v(\d+)';
[sim_todo, filecount] = all_topologies_simulate_count_todo(...
    max_trials, folder, pattern);
%}

%% loop over all phases
done = zeros(3,3,3,3); % keeps track of which topologies have been found already (up to symmetry)
M = [0 1 -1]; % index to interaction
count = 0;
for k=1:3^4
    disp(k);
    
    [i11, i12, i21, i22] = ind2sub([3, 3, 3, 3], k);
    
    % matrix associated with indices
    M_int = [M(i11) M(i12); M(i21) M(i22)];
    gM = [i22 i21; i12 i11];
    if done(i11,i12,i21,i22)
    	continue
    end
    count = count+1;
    done(i11,i12,i21,i22) = 1;
    done(gM(1,1),gM(1,2),gM(2,1),gM(2,2))=1;
    
    % ---- Generate random parameter set ----
    % bounds
    K_b = [1 10^3]; 
    Con_b = [1 10^3];
    
    % Latin hypercube
    idxK = find(M_int~=0);
    idxCon = find(sum(abs(M_int), 1)~=0);
    
    nK = numel(idxK); %sum(sum(abs(M_int))); 
    nCon = numel(idxCon); %sum(sum(abs(M_int), 1)>0);
    
    if nK+nCon==0
        continue
    end
    
    x = lhsdesign(n_pset, nK+nCon);
    
    % Visualize parameters
    %{
    figure;
    hold on
    xlim(K_b);
    ylim(Con_b);
    %}
    for idx1=1:n_pset
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
            try 
                load(fname);
            catch ME
                if strcmp(ME.identifier, 'MATLAB:load:couldNotReadFile')
                    warning('Could not find parameters.mat file for topology %d, pset %d',...
                        k, idx1);
                    break
                end
            end
        else
            save(fname, 'thisK', 'thisCon');
        end
        
        % count # simulations to do
        pattern = 'all_topologies_simulate-v(\d+)';
        [sim_todo, ~] = batch_sim_all_topologies_count_todo(...
            nsim, save_folder, pattern);
        
        % simulate trajectory
        %
        for count=1:sim_todo
            %disp(count);
            [cells_hist, period, t_onset] = time_evolution_save_func(N, a0,...
                Rcell, lambda, hill, noise, M_int, thisK, thisCon, Coff,...
                dist, InitiateI, p0, I0, tmax, save_folder);
        end
        %}
    end
end