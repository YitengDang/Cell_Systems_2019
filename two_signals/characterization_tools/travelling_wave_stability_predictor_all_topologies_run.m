% Predicts whether a travelling wave can propagate according to nearest
% neighbour interactions
% Apply computational search over all topologies to try to find parameters
% that allow for propagation of any kind of travelling wave
clear all
close all
set(0,'defaulttextinterpreter', 'latex')
%% Parameters
% Number of parameter sets to do
n_pset = 10^6;

% Manual input
gz = 15;
N = gz^2;
a0 = 1.5;
rcell = 0.2;
Rcell = rcell*a0;
lambda = [1 1.2];
%M_int = [0 1; -1 1]; % network 15 reversed
%{
M_int = [0 1; -1 1]; % network 15 reversed
M_int = [1 -1; 1 0]; % network 15
M_int = [1 1; -1 0]; % network 19
M_int = [1 -1; 1 1]; % network 33
M_int = [-1 -1; 1 1]; % network 34
M_int = [-1 1; -1 1]; % network 36
%}
%Con = [18 16];
%K = [0 9; 11 4];
%K = [4 11; 9 0];

% get pos, dist
mcsteps = 0;
[pos, dist] = initial_cells_random_markov_periodic(gz, mcsteps, rcell);

% Obtain from simulation
%{
folder = 'L:\BN\HY\Shared\Yiteng\two_signals\batch_sim_all_topologies_run2\selected';
subfolder = 'patterns\Network 33';
fname_str = 'tripple_wave_diagonally';
fname = fullfile(folder, subfolder, fname_str);
load(fname, 'save_consts_struct');

Con = save_consts_struct.Con;
K = save_consts_struct.K;
N = save_consts_struct.N;
gz = sqrt(N);
a0 = save_consts_struct.a0;
rcell = save_consts_struct.rcell;
Rcell = rcell*a0;
lambda = save_consts_struct.lambda;
M_int = save_consts_struct.M_int;

%}

% specify wave type and characteristics
wave_types_str = {'straight', 'inward bend', 'outward bend'};
wave_type = 1;
num_waves = 1; % number of waves
bandwidth = 1; % width of band of cells of one type

% order of states: F, M, B, E
%states_perm = [3 4 2 1]; % network 15 
%states_perm = [2 4 3 1]; % network 15 reversed
%{
states_perm = [2 4 3 1]; % network 15 reversed
states_perm = [3 4 2 1]; % network 15 
states_perm = [4 3 1 2]; % network 19
states_perm = [3 4 2 1]; % network 33
states_perm = [4 2 1 3]; % network 33/34
states_perm = [2 4 3 1]; % network 36
default_states = [0 0; 0 1; 1 0; 1 1];
%}

% save folder
%save_folder = 'L:\BN\HY\Shared\Yiteng\two_signals\trav_wave_stability_general\run2';
%save_folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\trav_wave_stability_general\run2';
save_folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\trav_wave_stability_general\run3_vary_a0_lambda12';
%% do for a fixed network & wave
%{
trav_wave_conds_met = travelling_wave_stability_predictor_general_func(gz, a0,...
    dist, rcell, lambda, M_int, Con, K,...
    wave_type, states_perm, num_waves, bandwidth);
%}
%% Loop over all waveforms
P = perms(1:4);
n_networks = 44;
for idx_P=12:size(P, 1)
    states_perm = P(idx_P, :);

    % loop over all phases
    done = zeros(3,3,3,3); % keeps track of which topologies have been found already (up to symmetry)
    M = [0 1 -1]; % index to interaction
    count = 0;
    networks_all = 1:3^4;
    networks_idx = [];
    
    trav_wave_cond_met = zeros(n_networks, n_pset);
    Con_all = zeros(n_networks, n_pset, 2);
    K_all = zeros(n_networks, n_pset, 2, 2);
    a0_all = zeros(n_networks, n_pset);
    lambda2_all = zeros(n_networks, n_pset);
    
    for k=networks_all
        fprintf('Waveform %d, Network %d \n', idx_P, k);
        [i11, i12, i21, i22] = ind2sub([3, 3, 3, 3], k);

        % matrix associated with indices
        M_int = [M(i11) M(i12); M(i21) M(i22)];
        gM = [i22 i21; i12 i11];
        if done(i11,i12,i21,i22)
            continue
        end
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

        % update network count
        count = count+1;
        networks_idx(count) = k;

        % Latin hypercube
        x = lhsdesign(n_pset, nK+nCon+2);
        
        % Visualize parameters
        %{
        figure;
        hold on
        xlim(K_b);
        ylim(Con_b);
        %}
        parfor idx1=1:n_pset
            thisK = zeros(2);
            thisCon = zeros(1,2);
            thisK(idxK) = (K_b(2) - K_b(1))*x(idx1, 1:nK) + K_b(1);
            thisCon(idxCon) = (Con_b(2) - Con_b(1))*x(idx1, nK+1:nK+nCon) + Con_b(1); 
            thisa0 = 10*x(idx1, nK+nCon+1);
            thislambda = [1 2*x(idx1, nK+nCon+2);];
            
            Con_all(count, idx1, :) = thisCon;
            K_all(count, idx1, :, :) = thisK;
            a0_all(count, idx1) = thisa0;
            lambda2_all(count, idx1) = thislambda(2);
            
            trav_wave_conds_met(count, idx1) = travelling_wave_stability_predictor_general_func(...
                gz, thisa0, dist, rcell, thislambda, M_int, thisCon, thisK,...
                wave_type, states_perm, num_waves, bandwidth);
            %}
        end
    end
    % fname
    fname_str = sprintf('trav_wave_conditions_check_wave_num_%d_type_%d_states_%d_%d_%d_%d',...
        num_waves, wave_type, states_perm(1), states_perm(2), states_perm(3), states_perm(4));
    fname = fullfile(save_folder, fname_str);
    save(fname, 'gz', 'a0', 'rcell', 'lambda', 'M_int',...
        'trav_wave_conds_met', 'Con_all', 'K_all', 'networks_idx');
end