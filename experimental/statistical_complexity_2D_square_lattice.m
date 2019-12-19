% Calculate statistical complexity of datasets of simulations of 2D square
% lattice systems (cellular automata)
clear all
close all
clc
%% Load data
% Single data set
load_folder = 'N:\tnw\BN\HY\Shared\Yiteng\cyclic_ca_2D\batch_2019_05_22';
fname_str = 'Cyclic_CA_l300_kappa4_r1_T1_sim_1';
load( fullfile(load_folder, fname_str) );
t_max = length(cells_hist)-1;
l = size(cells_hist{1}, 1);
L = l^2;
nbh_size = (2*r+1)^2; % size of neighborhood

% Multiple data sets

%% Replay trajectory
% Plot lattice
h=figure;
set(h, 'Position', [100 100 800 700]);
p=imagesc(cells_hist{1});
pbaspect([1 1 1])
colormap(parula(kappa))
c=colorbar;
title(sprintf('t=%d', 0));
set(c, 'YTick', []);
pause(0.1);

for t=1:t_max
    p.CData = cells_hist{t+1};
    title(sprintf('t=%d', t));
    pause(0.1);
end
%% Obtain light cones and conditional distributions
C_all = zeros(t_max-1, 1); % complexity vs time
for t=1:t_max-1
    disp(t);
    t_idx = t+1;
    
    cells_past = cells_hist{t_idx-1};
    cells_future = cells_hist{t_idx+1};
    
    past_cones = zeros(L, 1);
    future_cones = zeros(L, 1);
    
    %P_cond = zeros(kappa^nbh_size); % all states: not feasible
    
    for cell_i = 1:L
        % get cell neighborhood (linear indices)
        cells_nei = get_neighbors(cell_i, l, r);
        
        % Obtain past cones
        cells_nei_past = cells_past(cells_nei);
        this_past_cone = sum(cells_nei_past(:)'.*(2.^(0:nbh_size-1))); % conver to single number
        past_cones(cell_i) = this_past_cone;
        
        % Obtain future cones
        cells_nei_fut= cells_future(cells_nei);
        this_future_cone = sum(cells_nei_fut(:)'.*(2.^(0:nbh_size-1))); % conver to single number
        future_cones(cell_i) = this_future_cone;
        
        % Update conditional distribution     
        % (creates too large arrays, not feasible)
        %P_cond(this_future_cone, this_past_cone) = P_cond(this_future_cone, this_past_cone)+1;
    end
        
    %% Create a conditional distribution with only the observed
    % neighborhoods (cones)
    past_cones_uniq = unique(past_cones);
    future_cones_uniq = unique(future_cones);    
    P_cond = zeros(numel(past_cones_uniq), numel(future_cones_uniq));
    for idx = 1:L
        this_past_cone = past_cones(idx);
        this_future_cone = future_cones(idx);
        past_idx = find(past_cones_uniq == this_past_cone);
        future_idx = find(future_cones_uniq == this_future_cone);
        P_cond(past_idx, future_idx) = P_cond(past_idx, future_idx)+1;
    end
    % normalize P_cond
    P_cond = P_cond./sum(P_cond, 2); 

    %% Construct causal states
    P_cond_uniq = {}; % unique 
    P_causal_states = []; % probability distribution of causal states
    % shuffle states (to do)
    for idx=1:numel(past_cones_uniq)
        this_P_cond = P_cond(idx, :);

        % check for each distribution whether it matches the current
        match = 0;
        match_idx = 1;
        for idx2=1:numel(P_cond_uniq)
            that_P_cond = P_cond_uniq{idx2};
            if all(this_P_cond==that_P_cond)
                match = 1;
                match_idx = idx2;
                break
            end
        end

        if match
            % add to causal state count
            P_causal_states(idx2) = P_causal_states(idx2) + 1;
        else 
            % create a new causal state
            P_cond_uniq{end+1} = this_P_cond;
            P_causal_states(end+1) = 1;
        end
    end

    % Normalize causal state distribution
    P_causal_states = P_causal_states./sum(P_causal_states);

    %% Calculate complexity
    C_all(t) = -sum( P_causal_states .*log2(P_causal_states));
end

%% Save analyzed trajectory
save_folder = 'N:\tnw\BN\HY\Shared\Yiteng\cyclic_ca_2D\batch_2019_05_22_analyzed';
fname_str_out = strcat(fname_str, '_analyzed_complexity');
save(fullfile(save_folder, fname_str_out), 'C_all');

%% Plot complexity vs time

figure;
plot(1:t_max, C_all, 'o-', 'LineWidth', 2);
xlabel('Time');
ylabel('Stat. complexity (bits per site)');

%% 
function cells_nei = get_neighbors(cell_i, l, r)
    [dx, dy] = meshgrid(-r:r, -r:r); % Moore neighborhood
    [cell_x, cell_y] = ind2sub([l l], cell_i);
    x_nei = mod(cell_x+dx-1,l)+1; % x-indices of neighbours
    y_nei = mod(cell_y+dy-1,l)+1; % y-indices of neighbours
    cells_nei = sub2ind([l l], x_nei, y_nei); % linear indices of cells
end