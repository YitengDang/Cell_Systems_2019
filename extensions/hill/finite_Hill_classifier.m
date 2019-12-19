%% Classify finite Hill trajectories by final state
clear all
close all

%% Test single trajectory
load_folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\trav_wave_hill\TW_formation_network_15\Parameter_set_1';
fname_str = 'two_signal_mult_N225_params_1_sim_1_hill_10_vpa_digits_36_t_out_10000_period_Inf_vpa_digits_36-v1';
load(fullfile(load_folder, fname_str), 'cells_hist', 'period', 't_out');
N = size(cells_hist{1}, 1);
t_max = 10^4;

[static, homogeneous] = trajectory_classifier(cells_hist, period);
t_max_reached = (t_out==t_max);
fprintf("Static = %d, homogeneous = %d, t_max reached = %d \n",...
    static, homogeneous, t_max_reached);

% move selected files
%{
class_id = '';
if static && ~homogeneous
    class_id = 'static_pattern';
elseif ~static && homogeneous
    class_id = 'collective_oscillation';
elseif ~static && ~homogeneous
    class_id = 'dynamic_spatial_pattern';
end

if numel(class_id)>0
    disp('Copying file...');
    old_file = fullfile(load_folder, strcat(fname_str, '.mat'));
    new_file = fullfile(parent_folder, 'simulations_sorted', class_id, strcat(fname_str, '.mat'));
    [status, msg] = copyfile(old_file, new_file);
    if status
        fprintf('File copied to %s \n', new_file);
    else
        error(msg);
    end
end
%}
%% Classify batch simulations, move files
% settings
gz = 15;
N = gz^2;
vpa_digits = 36;
t_max = 10^4;
num_params = 10;
nruns = 100;
hill_all = [1 2 5 10 Inf];
parent_folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\trav_wave_hill\TW_formation_network_15';

% variables to store
static_all = zeros(num_params, nruns, numel(hill_all));
homogeneous_all = zeros(num_params, nruns, numel(hill_all));
filecount_all = zeros(num_params, nruns, numel(hill_all));
t_max_reached_all = zeros(num_params, nruns, numel(hill_all));

%% loop over parameters
for idx_param=2:num_params
    %idx_param = 1;

    % get all files
    subfolder = sprintf('Parameter_set_%d', idx_param);
    folder = fullfile(parent_folder, subfolder);
    listing = dir(folder);
    num_files = numel(listing)-2; %first two entries are not useful
    count = 0;
    for i = 1:num_files
        filename = listing(i+2).name;
        % remove extension and do not include txt files
        [~,name,ext] = fileparts(filename);
        if strcmp(ext, '.mat')
            count = count + 1;
            names{count} = name;
        end
    end

    % pattern
    pattern = sprintf('two_signal_mult_N%d_params_%d_sim_%s_hill_%s_vpa_digits_%d_t_out_%s_period_%s_vpa_digits_%d-v%s',...
        N, idx_param, '(\d+)', '(\d+)', vpa_digits, '(\d+)', '(\d+|Inf)', vpa_digits, '(\d+)');
    % classify all files
    for i=1:numel(names)
        [tokens, ~] = regexp(names{i}, pattern, 'tokens', 'match');
        if ~isempty(tokens)
            disp(names{i});
            idx_sim = str2double(tokens{1}{1});
            hill = str2double(tokens{1}{2});
            idx_hill = find(hill == hill_all, 1);

            % classify final state
            load(fullfile(folder, names{i}), 'cells_hist', 'period', 't_out');
            [static, homogeneous] = trajectory_classifier(cells_hist, period);
            t_max_reached = (t_out==t_max);

            % update variables
            filecount_all(idx_param, idx_sim, idx_hill) = filecount_all(idx_param, idx_sim, idx_hill)+1;
            static_all(idx_param, idx_sim, idx_hill) = static;
            homogeneous_all(idx_param, idx_sim, idx_hill) = homogeneous;
            t_max_reached_all(idx_param, idx_sim, idx_hill) = t_max_reached;

            % move selected files (optional)
            %{
            class_id = '';
            if static && ~homogeneous
                class_id = 'static_pattern';
            elseif ~static && homogeneous
                class_id = 'collective_oscillation';
            elseif ~static && ~homogeneous
                class_id = 'dynamic_spatial_pattern';
            end
            if ~isempty(class_id) && ~t_max_reached
                old_file = fullfile(parent_folder, subfolder, strcat(names{i}, '.mat'));
                new_file = fullfile(parent_folder, 'simulations_sorted', class_id, strcat(names{i}, '.mat'));
                [status, msg] = copyfile(old_file, new_file);
                if status
                    fprintf('File copied to %s \n', new_file);
                else
                    error(msg);
                end
            end
            %}
        end
    end

end

% Save analyzed data
save_folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\trav_wave_hill\TW_formation_network_15';
fname_str = sprintf('analyzed_data_classification_vs_hill_num_params_%d_nruns_%d_digits_%d',...
    num_params, nruns, vpa_digits);
save(fullfile(save_folder, fname_str), 'filecount_all', 'static_all', 'homogeneous_all');

%% Also load infinite Hill data
% -- to do--

%% Plot analyzed data and save
% classification: 1=dynamic spatial, 2=collective oscillation, 3=static pattern
% 4=no pattern, 5=max. simulation time reached
save_path_fig = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\trav_wave_hill\TW_formation_network_15';

classification_all = filecount_all + static_all.*2 + homogeneous_all;
classification_all(t_max_reached_all==1) = 5;

% plot for all parameter sets separately
frac_all = zeros(num_params, numel(hill_all)-1, 5); % five classes
h=figure;
hold on
for idx_param=1:num_params
    % idx_param = 1;
    cats = sprintfc('%d', 1:5);
    for idx_hill=1:numel(hill_all)-1
        data = squeeze(classification_all(idx_param, :, idx_hill));
        [X, ~] = histcounts( categorical(data), cats );
        frac_all(idx_param, idx_hill, :) = X/(nruns);
    end
    
    subplot(3,4,idx_param)
    bar( squeeze(frac_all(idx_param,:,:)), 'stacked' )
    title(sprintf('Pset %d', idx_param));
    %legend(cats);
    xlabel('Hill coefficient');
    ylabel('Fraction')
    set(gca, 'FontSize', 20);
    
    xticklabels = sprintfc('%d', hill_all(1:end-1) );
    set(gca, 'xticklabels', xticklabels);
end

% Save figure
qsave = 0;
fname = fullfile(save_path_fig, strcat('analyzed_data_', subfolder,...
    sprintf('_nruns_%d_digits_%d', nruns, digits), '_classification_separate_psets'));
save_figure(h, 12, 8, fname, '.pdf', qsave);
%% Stack results from different psets together
frac_all_merged = squeeze(sum(frac_all, 1))./num_params;
h=figure;
bar( frac_all_merged, 'stacked' )
xlabel('Hill coefficient');
ylabel('Fraction')
set(gca, 'FontSize', 32);
legend({'dynamic spatial pattern', 'collective oscillation', 'static pattern',...
    'no pattern', 't_{max} reached'}, 'FontSize', 20);

xticklabels = sprintfc('%d', hill_all(1:end-1) );
set(gca, 'xticklabels', xticklabels);

% Save figure
qsave = 0;
fname = fullfile(save_path_fig, strcat('analyzed_data_', subfolder,...
    sprintf('_nruns_%d_digits_%d', nruns, digits), '_classification_merged_psets'));
save_figure(h, 12, 8, fname, '.pdf', qsave);

%% Analyse simulation times

%% Functions
function [static, homogeneous] = trajectory_classifier(cells_hist, period)
    % (I) static or dynamic?
    if period<Inf
        % check whether the periodicity corresponds to a "real" oscillation
        % by checking for each cell whether it "oscillates" or stays constant
        cells_one_period = cell2mat(cells_hist(end-period:end)); % size: N * (period+1)*2;

        eps_I = 10^(-5);
        cond1 = max(range(cells_one_period(:, 1:2:end-1), 2))<eps_I; % range of values of gene 1
        cond2 = max(range(cells_one_period(:, 2:2:end), 2))<eps_I; % gene 2
        static = (cond1 && cond2);
    else
        cells_one_period = cells_hist{end};
        static = 1;
    end
    % (II) homogeneous or spatially structured?
    if static
        % static patterns: check final state only
        cells_final = cells_one_period(:, end-1:end); %use previous variable, select only last time point
    else
        % dynamic patterns: check all states of one period
        cells_final = cells_one_period; %use previous variable, select only last time point
    end

    eps_II = 10^(-5);
    homogeneous = all(var(cells_final, 0, 1) < eps_II);
end