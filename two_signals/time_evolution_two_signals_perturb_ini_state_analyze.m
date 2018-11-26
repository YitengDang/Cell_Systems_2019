clear all
close all
%%
% Parameters
N = 225;
K12 = 8;
num_cells_changed = 10;

% Load data

% filter on original data first
pattern = sprintf('^two_signal_mult_N%d_K12_%d_t_out_%s_period_%s-v%s$',...
	N, K12, '(\d+)', '(\d+|Inf)', '(\d+)');

parent_folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\sweep K12 new lattice\sensitivity_init_cond';
subfolder = sprintf('K12_%d', K12);
listing = dir(fullfile(parent_folder, subfolder));
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

%% then, consider the data with cells flipped
names_flipped = cell(numel(names_orig), 1);

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
        [~,fname,ext] = fileparts(filename);
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
end

%% Load data and process it
trajectory_distances = cell(size(names_flipped));
for i=1:numel(names_orig)
    % initialize empty cell array
    trajectory_distances{i} = {};
    
    % load original data
    load(fullfile(parent_folder, subfolder, names_orig{i}), 'cells_hist');
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
%% Plot distance trajectories
for i=1:numel(trajectory_distances)
    h = figure;
    hold on
    for i2=1:numel(trajectory_distances{i})
        cells_dist_temp = trajectory_distances{i}{i2};
        t_data = 0:numel(cells_dist_temp)-1;
        plot(t_data, cells_dist_temp);
    end
    ylim([0 1]);
    xlabel('time');
    ylabel('distance');
    set(gca, 'FontSize', 20);
    
    % save figure
    qsave = 1;
    save_folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\sweep K12 new lattice\sensitivity_init_cond\figures';
    fname_str = sprintf('%s_%d_cells_changed_distance_vs_t_plot', names_orig{i},  num_cells_changed);
    fname = fullfile(save_folder, fname_str);
    save_figure(h, 10, 8, fname, '.pdf', qsave);
end