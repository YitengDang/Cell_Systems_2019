% Calculate statistical complexity of datasets
clear all
close all
clc
%% Load data

% Single data set
load_folder = 'H:\My Documents\Multicellular automaton\app\data\time_evolution\test_data_statistical_complexity';
fname_str = 'uniform_lattice_osc_period_2';
load( fullfile(load_folder, fname_str) );
t_max = length(cells_hist)-1;
N = size(cells_hist{1}, 1);

% Multiple data sets

%% Obtain conditional distributions
P_cond_all = cell(N, 1); % 1=(0,0), 2=(0,1), 3=(1,0), 4=(1,1)
past_states_all = cell(N, 1);
for cell_i = 1:N
    past_states = []; % past states l- that have occured so far
    past_state_counts = []; % counts of the past states l- in "past_states"
    P_cond_cell = {}; % all P(l+|l-) for cell_i
    % Calculate transition probabilities from data at different time points
    % (ergodic assumption)
    for t=1:t_max-1
        t_idx = t+1;
        cells_past = cells_hist{t_idx-1};
        cells_future = cells_hist{t_idx+1};
        
        % Obtain past cone
        % -- code to get all cell neighbors
        cells_past_cone = [];
        %cell_idx_past = cells_past(cell_i, :)*[2; 1];
        %states_count(cell_idx_past) = states_count(cell_idx_past) + 1;
        this_l_min = cells_past*(2.^(0:6));
        
        % Obtain future cone
        %cell_idx_future = cells_future(cell_i, :)*[2; 1];
        this_l_plus = cells_future*(2.^(0:6));
        
        % check if past state has occured before
        % convert past_cone_states to number
        idx_temp = (past_states == this_l_min);
        if ~isempty(idx_temp) % Existing past state
            past_state_counts(idx_temp) = past_state_counts(idx_temp) + 1;
            % Update P_cond
            this_P_cond = P_cond_cell{idx_temp};
            this_P_cond(this_l_plus, this_l_min) = this_P_cond(this_l_plus, this_l_min) + 1;
            P_cond_cell{idx_temp} = this_P_cond;
        else % New past state observed
            past_state_end(end+1) = this_l_min;
            past_state_counts(end+1) = 1;
            % Add new P_cond to list
            this_P_cond = zeros(2^7);
            this_P_cond(this_l_plus, this_l_plus) = 1;
            P_cond_cell{end+1} = this_P_cond;
            % P_cond_cell{end+1} = zeros(2^7);
        end
    end
    % normalize P_cond
    sel_idx = states_count>0; % only select indices of states that are present
    P_cond(sel_idx, :) = P_cond(sel_idx,:)./states_count(sel_idx); 
    
    if ~(all(sum(P_cond, 2)==0 | sum(P_cond, 2)==1))
        error('Wrong probabilities')
    end
    
    % store results
    P_cond_all(cell_i,:,:) = P_cond;
    states_count_all(cell_i, :) = states_count;
end

%% Construct causal states

% 
cell_i = 1;
P_cond = squeeze(P_cond_all(cell_i,:,:));

%% Calculate complexity