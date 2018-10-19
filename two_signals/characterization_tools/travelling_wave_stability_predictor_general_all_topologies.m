% Predicts whether a travelling wave can propagate according to nearest
% neighbour interactions
% Apply computational search over all topologies to try to find parameters
% that allow for propagation of any kind of travelling wave
clear all
close all
set(0,'defaulttextinterpreter', 'latex')
%% Parameters
% Number of parameter sets to do
n_pset = 10000;

% Manual input
gz = 15;
N = gz^2;
a0 = 1.5;
rcell = 0.2;
Rcell = rcell*a0;
lambda = [1 1.2];
M_int = [1 1; -1 0]; % network 19 
%{
M_int = [0 1; -1 1]; % network 15 reversed
M_int = [1 -1; 1 0]; % network 15
M_int = [1 1; -1 0]; % network 19
M_int = [1 -1; 1 1]; % network 33
M_int = [-1 -1; 1 1]; % network 34
M_int = [-1 1; -1 1]; % network 36
%}
%Con = [18 16];
%K = zeros(2);
%K = [0 9; 11 4];

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
states_perm = [4 2 1 3]; 
%{
states_perm = [2 4 3 1]; % network 15
states_perm = [4 3 1 2]; % network 19
states_perm = [3 4 2 1]; % network 33
states_perm = [4 2 1 3]; % network 33/34
states_perm = [2 4 3 1]; % network 36
%}

% fname
fname_str = sprintf('');
%% loop over all phases
done = zeros(3,3,3,3); % keeps track of which topologies have been found already (up to symmetry)
M = [0 1 -1]; % index to interaction
count = 0;
networks_all = 1:3^4;
trav_wave_cond_met = zeros(numel(networks_all), n_pset);

for k=networks_all
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
        
        % calculate stability
        
        % --> wrap in function
        trav_wave_cond_met(l, idx1) = 1; % output
        %}
    end
end

