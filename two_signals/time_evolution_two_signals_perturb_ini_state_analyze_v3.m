%% Analyze the results from perturbing cell states
clear all
close all
clc
%%
% Parameters
N = 225;
K12 = 9;
num_cells_changed = 1;

%% Load data
%
% filter on original data first
pattern = sprintf('^two_signal_mult_N%d_K12_%d_t_out_%s_period_%s-v%s$',...
	N, K12, '(\d+)', '(\d+|Inf)', '(\d+)');
parent_folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\sweep K12 new lattice\sensitivity_init_cond';
subfolder = sprintf('K12_%d_%dcell_changed', K12, num_cells_changed);
subfolder2 = 'originals';
listing = dir(fullfile(parent_folder, subfolder, subfolder2));
num_files = numel(listing)-2; %first two entries are not useful
count = 0;
names_orig = {};
for i = 1:num_files
    filename = listing(i+2).name;
    % remove extension and do not include txt files
    [~,fname,ext] = fileparts(filename);
    if strcmp(ext, '.mat')
        [~, tokens] = regexp(fname, pattern, 'match', 'tokens');
        if ~isempty(tokens)
            disp(fname);
            count = count + 1;
            names_orig{count} = fname;
        end
    end
end

%% Then, consider the data with cells flipped
names_flipped = cell(numel(names_orig), 1);

% get list of filenames
listing = dir(fullfile(parent_folder, subfolder));
num_files = numel(listing)-2;
% get and organize filenames of new data
for i=1:numel(names_orig)
    name_orig = names_orig{i};
    names_flipped{i} = {};
    fprintf('old file = %s \n', name_orig);
    
    %pattern2 = sprintf('two_signal_mult_N%d_K12_%d_t_out_%d_period_%d-v%d_%dcells_changed-v%s',...
    %	N, K12, '(\d+)', '(\d+|Inf)', '(\d+)', num_cells_changed, '(\d+)');
    pattern2 = sprintf('%s_%dcells_changed-v%s', name_orig, num_cells_changed, '(\d+)');
    
    count = 0;
    for i2=1:num_files
        %disp(i2);
        filename = listing(i2+2).name;
        % remove extension and do not include txt files
        [~, fname, ext] = fileparts(filename);
        %disp(fname);
        if strcmp(ext, '.mat')
            [~, tokens] = regexp(fname, pattern2, 'match', 'tokens');
            if ~isempty(tokens)
                % correct data format -> load file, store results
                %disp(fname);
                count = count+1;
                names_flipped{i}{count} = fname;
            end
        end
    end
    fprintf('Count = %d \n', count);
end

%% Get number of names
sizes = zeros(numel(names_flipped), 1);
for i=1:numel(names_flipped)
    sizes(i) = size(names_flipped{i}, 2);
end
if numel(unique(sizes))>1
    error('Simulations per file not equal!');
else
    nsim = sizes(1);
end
%% Load data and process it
trajectory_distances = cell(size(names_flipped));
for i=1:numel(names_orig)
    % initialize empty cell array
    trajectory_distances{i} = {};
    
    % load original data
    load(fullfile(parent_folder, subfolder, subfolder2, names_orig{i}), 'cells_hist');
    cells_hist_orig = cells_hist;
    
    % look at flipped trajectories
    for i2=1:numel(names_flipped{i})
        fname = names_flipped{i}{i2};
        disp(fname);
        
        % load data
        load(fullfile(parent_folder, subfolder, fname), 'cells_hist');
        
        % find distances over all times
        t_final = min( numel(cells_hist), numel(cells_hist_orig) );
        cells_dist_temp = zeros(t_final, 1);
        for t=1:t_final
            cells_diff = ~all(cells_hist_orig{t}==cells_hist{t}, 2); %not both genes the same
            cells_dist_temp(t) = sum(cells_diff)/N;
        end
        trajectory_distances{i}{i2} = cells_dist_temp;
    end
end
%% Save analyzed data
save_folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\sweep K12 new lattice\sensitivity_init_cond';
fname_str = sprintf('analyzed_data_K12_%d_%druns', K12, nsim);
save( fullfile(save_folder, fname_str), 'names_orig', 'names_flipped', 'trajectory_distances'); 
%}
%% Load analyzed data
save_folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\sweep K12 new lattice\sensitivity_init_cond';
fname_str = sprintf('analyzed_data_K12_%d_%druns', K12, nsim);
load( fullfile(save_folder, fname_str), 'names_orig', 'names_flipped', 'trajectory_distances'); 

%%
%{
% load original data
load(fullfile(parent_folder, subfolder, names_orig{i}), 'cells_hist');
cells_hist_orig = cells_hist;
    
i = 1;
i2 = 1;
fname = names_flipped{i}{i2};
disp(fname);

% load data
load(fullfile(parent_folder, subfolder, fname), 'cells_hist');

% find distances over all times
t_final = min( numel(cells_hist), numel(cells_hist_orig) );
cells_dist_temp = zeros(t_final, 1);
for t=1:t_final
    cells_diff = ~all(cells_hist_orig{t}==cells_hist{t}, 2); %not both genes the same
    cells_dist_temp(t) = sum(cells_diff)/N;
end
%}

%% Get average distance
avg_distances = cell( numel(trajectory_distances), 1 );
for i=1:numel(trajectory_distances)
    % get longest simulation time
    tmax = 0;
    for i2=1:numel(trajectory_distances{i})
        tmax = max( tmax, numel(trajectory_distances{i}{i2})-1 );
    end
    
    % get average distance
    avg_distances_temp = zeros(tmax+1, 1); % avg distance
    avg_dist_traj_count = zeros(tmax+1, 1); % count # trajectories
    for i2=1:numel(trajectory_distances{i})
        cells_dist_temp = zeros(tmax+1, 1);
        t_temp = trajectory_distances{i}{i2};
        cells_dist_temp(1:numel(t_temp)) = t_temp;
        
        %if tmax+1-numel(cells_dist_temp) >0
        %    cells_dist_temp = padarray(cells_dist_temp, tmax-numel(cells_dist_temp));
        %end
        avg_distances_temp = avg_distances_temp + cells_dist_temp;
        avg_dist_traj_count(1:numel(t_temp)) = avg_dist_traj_count(1:numel(t_temp)) + 1;
    end
    avg_distances_temp = avg_distances_temp./avg_dist_traj_count; %numel(trajectory_distances{i});
    
    % store result
    avg_distances{i} = avg_distances_temp;
end

%% Plot distance trajectories
%
for i=1:numel(trajectory_distances)
    figure;
    this_tmax = min(tmax, numel(avg_distances{i})-1);
    plot(0:this_tmax, avg_distances{i}(1:this_tmax+1), 'r', 'LineWidth', 2);
    hold on
    for i2=1:numel(trajectory_distances{i})
        cells_dist_temp = trajectory_distances{i}{i2};
        t_data = 0:numel(cells_dist_temp)-1;
        plot(t_data, cells_dist_temp, 'Color', [0.6 0.6 0.6 0.3]);
    end
    ylim([0 1]);
    xlabel('time');
    ylabel('distance');
    set(gca, 'FontSize', 20);
    
    % save figure
    qsave = 1;
    save_folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\sweep K12 new lattice\sensitivity_init_cond\figures';
    fname_str = sprintf('%s_%d_cells_changed_distance_vs_t_plot_w_avg_%d_sims',...
        names_orig{i}, num_cells_changed, nsims);
    fname = fullfile(save_folder, fname_str);
    save_figure(h, 10, 8, fname, '.pdf', qsave);
end
%}
%% Plot specified range 
%
%i = 5;
tmax = 60;
for i=1:10
    h = figure;
    plot(0:tmax, avg_distances{i}(1:tmax+1), 'r', 'LineWidth', 2);

    hold on
    for i2=1:numel(trajectory_distances{i})
        cells_dist_temp = trajectory_distances{i}{i2};
        %t_data = 0:numel(cells_dist_temp)-1;
        t_range = 0:min(tmax, numel(cells_dist_temp)-1);
        %plot(t_range, cells_dist_temp(t_range+1) );
        plot(t_range, cells_dist_temp(t_range+1), 'Color', [0.6 0.6 0.6 0.3]);
    end
    ylim([0 1]);
    xlabel('time');
    ylabel('distance');
    set(gca, 'FontSize', 20);
    
    % save figure
    qsave = 1;
    save_folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\sweep K12 new lattice\sensitivity_init_cond\figures';
    %fname_str = sprintf('%s_%d_cells_changed_distance_vs_t_plot_t_range_%dto%d',...
    %    names_orig{i},  num_cells_changed, t_range(1), t_range(end) );
    fname_str = sprintf('%s_%d_cells_changed_distance_vs_t_plot_t_range_%dto%d_w_avg_%d_sims',...
        names_orig{i}, num_cells_changed, t_range(1), t_range(end), nsims );
    fname = fullfile(save_folder, fname_str);
    save_figure(h, 10, 8, fname, '.pdf', qsave);
end
%}