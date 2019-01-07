% Analyze batch simulation results for moving cells
close all
clear all
% maxNumCompThreads(4);
% warning off
set(0, 'defaulttextinterpreter', 'latex');

%% Set folders
%parent_folder = 'H:\My Documents\Multicellular automaton\temp';
parent_folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\moving_growing_cells';

% Save data folder
save_folder = fullfile(parent_folder, 'figures');

% Save figures folder
save_fig_folder = fullfile(parent_folder, 'figures');
%% Load simulations for moving cells
gz = 16;
N = gz^2;
rcell_sigma = 0.1;
k_growth = 1.5;
sigma_D = 0.1;
t_out = 1000;

%subfolder = 'Network_15';
subfolder = 'sustained_inhomogeneity';
folder = fullfile(parent_folder, subfolder);

% Load growing/moving cells simulations
nruns = 50;
p_final_all = zeros(nruns, 2);
period_all = zeros(nruns, 1);
% plot p vs t all together
%{
h1 = figure;
hold on
h2 = figure;
hold on
%}
for jj=1:nruns
    %fname_str = 'one_signal_sigma_D_0p001_t_out_1000-v1';
    fname_str = strrep(sprintf('two_signals_rcell_sigma_%.1f_K_growth_%.1f_sigma_D_%.3f_t_out_%d-v%d',...
        rcell_sigma, k_growth, sigma_D, t_out, jj), '.' ,'p');
    fname = fullfile(folder, strcat(fname_str, '.mat'));
    if exist(fname, 'file')==2
        disp(fname);
        load(fname, 'save_consts_struct', 'cells_hist', 't_out',...
            'changed', 'positions', 'distances', 'positions_all', 'rcell_hist');
        
        % Store final p
        cells_final = cells_hist{end};
        p_final_all(jj, :) = sum(cells_final)/N;      
        
        % plot p vs t
        pt = zeros( numel(cells_hist), 2 );
        for tt=1:numel(cells_hist)
            pt(tt, :) = mean(cells_hist{tt}, 1);
        end
        figure;
        hold on
        plot(0:numel(cells_hist)-1, pt(:,1), 'b--');
        plot(0:numel(cells_hist)-1, pt(:,2), 'r--');        
        
        % Determine period
        period = Inf;
        [period_ub, ~] = periodicity_test_short(cells_hist); 
        
        % if periodicity found, refine check to find period
        if period_ub<Inf
            t_check_init = 0;
            decimals = Inf;
            [period, t_onset] = periodicity_test_detailed(cells_hist, t_check_init,...
                period_ub, decimals);
            t_out = t_onset + period; 
        end
        period_all(jj) = period;        
    end
end

%% Load negative control simulations
count_by_type = zeros(3, 1); 
% type 1: reaches steady state, type 2: periodic trajectories, type 3:
% simulations reaching tmax / indeterminate

files = dir(folder);
files = files(3:end);
rcell_sigma = 0;
sigma_D = 0;
k_growth = 1.5;
tmax = 1000;

p_final_all_nc = [];

count = 0;
for ii=1:numel(files)
    filename = files(ii).name;
    %pattern = strrep(sprintf(...
    %    'two_signals_rcell_sigma_%.1f_K_growth_%.1f_sigma_D_%.3f_t_out_%s_neg_control-v%s',...
    %    rcell_sigma, k_growth, sigma_D, '(\d+)', '(\d+)'), '.', 'p');
    pattern = strrep(sprintf(...
        'two_signals_rcell_sigma_%.1f_sigma_D_%.3f_t_out_%s_period_%s_neg_control-v%s',...
        rcell_sigma, sigma_D, '(\d+)', '(\d+|Inf)', '(\d+)'), '.', 'p');
    [tokens, ~] = regexp(filename, pattern, 'tokens', 'match');
    if ~isempty(tokens)
        disp(filename);
        load(fullfile(folder, filename), 'save_consts_struct', 'cells_hist', 't_out',...
            'period', 'changed', 'positions', 'distances', 'positions_all', 'rcell_hist');
        count = count + 1;
        
        % classify trajectory
        if t_out==tmax % simulations reaching tmax
            count_by_type(3) = count_by_type(3)+1;
        elseif period==Inf % simulations reaching steady state
            count_by_type(1) = count_by_type(1)+1;
        else % t_out<tmax, period<Inf -> periodic trajectories
            count_by_type(2) = count_by_type(2)+1;
        end
        
        % store final p
        p_final_all_nc(count, :) = mean(cells_hist{end});
        
        % check for periodicity
        %{
        period = Inf;
        [period_ub, ~] = periodicity_test_short(cells_hist); 
        
        % if periodicity found, refine check to find period
        if period_ub<Inf
            t_check_init = 0;
            decimals = Inf;
            [period, t_onset] = periodicity_test_detailed(cells_hist, t_check_init,...
                period_ub, decimals);
            t_out = t_onset + period; 
        end
        
        % re-save as new file
        save_consts_struct.k_growth = 0;
        save_consts_struct.K_growth = 0;
                
        new_filename = strrep(sprintf('two_signals_rcell_sigma_%.1f_sigma_D_%.3f_t_out_%s_period_%d_neg_control-v%s',...
            rcell_sigma, sigma_D, tokens{1}{1}, period, tokens{1}{2}), '.', 'p');
        disp(new_filename);
        
        save(fullfile(folder, new_filename), 'save_consts_struct', 'cells_hist', 't_out', 'period',...
            'changed', 'positions', 'distances', 'positions_all', 'rcell_hist');
        %}
    end
    %}
end
%% Save analyzed data
%
save_fname_str = strrep(sprintf('analyzed_data_%s_t_max_%d_nruns_%d_rcell_sigma_%.1f_K_growth_%.1f_sigma_D_%.3f',...
    subfolder, tmax, nruns, rcell_sigma, k_growth, sigma_D), '.', 'p');
save_file = fullfile(parent_folder, save_fname_str);
save(save_file, 'subfolder', 'p_final_all', 'p_final_all_nc', 'period_all',...
    'rcell_sigma', 'k_growth', 'sigma_D', 't_out');
%}

%% Plot final p
% For moving cells
bins = 0:0.1:1;
h = figure;
hold on
histogram( p_final_all(:,1), bins );
histogram( p_final_all(:,2), bins );
title(sprintf('n=%d', nruns));
xlabel('$p_{final}$');
ylabel('Count');
set(gca, 'FontSIze', 20);

% For negative control
bins = 0:0.1:1;
h = figure;
hold on
histogram( p_final_all_nc(:,1), bins );
histogram( p_final_all_nc(:,2), bins );
title(sprintf('n=%d', nruns));
xlabel('$p_{final}$');
ylabel('Count');
set(gca, 'FontSIze', 20);

%% Plot classification
h = figure;
p=bar([count_by_type'; [0 0 nruns] ], 'stacked');
legend({'Steady state', 'Periodic', 'Unfinished'});
set(gca, 'XTickLabel', {'Stationary cells', 'Moving & growing cells'});
%xtickangle(45);
xlim([0.25 2.75]);
xlabel('');
ylabel('Count');
set(gca, 'FontSize', 20);

qsave = 1;
fname_str = 'Sustained_inhomogeneity_analyzed_data_bar';
path_out = fullfile(parent_folder, fname_str);
save_figure(h, 10, 8, path_out, '.pdf', qsave)