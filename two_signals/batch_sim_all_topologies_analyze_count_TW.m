%% Obtain statistics on all possible topologies by simulation
% Count fraction of travelling wave trajectories
% Based on pre-selection of trajectories
clear variables
close all
clc
set(0, 'defaulttextinterpreter', 'tex');
%set(0, 'defaulttextinterpreter', 'none');
%% Parameters and settings
% Settings
remote = 0;

% Note: increasing nsim at n_pset is always possible. However, increasing
% n_pset leads to data sets that do not form a perfect LHS sample
n_pset = 10000; % number of parameter sets to do
nsim = 10; % number of simulations per parameter set
tmax = 10000;

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

digits = -floor(log10(1/N)); % number of significant digits required to determine constancy of p

% get pos, dist
mcsteps = 0;
[pos, dist] = initial_cells_random_markov_periodic(gz, mcsteps, rcell);

% folders
load_folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\batch_sim_all_topologies_run2\selected';
save_folder = 'H:\My Documents\Multicellular automaton\figures\two_signals\batch_sim_all_topologies_run2\count_TW';
if remote
    load_folder = strrep(load_folder, 'N:\', 'W:\staff-bulk\');
    save_folder = strrep(save_folder, 'H:\', 'W:\staff-homes\d\yitengdang\');
end
%% Load full data
%{
networks_sel = [15 19 33 34 36]; %[15  16	19	20	32	33	34	36	43];
TW_count_strict = zeros( numel(networks_sel), n_pset );
TW_count_loose = zeros( numel(networks_sel), n_pset );
K_all_by_network = cell(numel(networks_sel), 1);
Con_all_by_network = cell(numel(networks_sel), 1);

for network_idx=1:numel(networks_sel) %networks_sel
    network = networks_sel(network_idx);
    subfolder = sprintf('Network_%d', network);
    
    K_all_temp = [];
    Con_all_temp = [];
    
    % Get filenames
    listing = dir( fullfile(load_folder, subfolder) );
    num_files = numel(listing)-2; %first two entries are not useful
    fprintf('Network %d, num_files = %d \n', network, num_files);
    
    for ii=1:num_files
        fprintf('File %d out of %d \n', ii, num_files);
        pattern = sprintf('selected_network_%d_param_%s_sim_%s_t_out_%s_period_%s',...
            network, '(\d+)', '(\d+)', '(\d+)', '(Inf|\d+)');
        filename = listing(ii+2).name;
        % remove extension and do not include other files
        [~,name,ext] = fileparts(filename);
        if strcmp(ext, '.mat')
            [~, tokens] = regexp(name, pattern, 'match', 'tokens');
            if ~isempty(tokens)
                % identify parameter set
                pset_idx = str2double(tokens{1}{1});
                
                % Select ones with correct period
                this_period = str2double(tokens{1}{4});
                if mod(this_period, gz)==0
                    % Check whether selected files are TWs
                    %disp('TW found!');
                    %disp(filename);
                    load(fullfile(load_folder, subfolder, filename),...
                        'cells_hist', 'save_consts_struct', 'distances');
                    a0 = save_consts_struct.a0;
                    [trav_wave, trav_wave_2] = travelling_wave_test(cells_hist, a0,...
                        this_period, numel(cells_hist)-1, distances, digits);
                    TW_count_strict(network_idx, pset_idx) = ...
                        TW_count_strict(network_idx, pset_idx) + trav_wave;
                    TW_count_loose(network_idx, pset_idx) = ...
                        TW_count_loose(network_idx, pset_idx) + trav_wave_2;
                    
                    % Store K, Con values
                    if trav_wave_2
                        count = sum(sum(TW_count_loose(network_idx, :, :)));
                    	K_all_temp(count, :, :) = save_consts_struct.K;
                        Con_all_temp(count, :) = save_consts_struct.Con;
                    end
                end
            end
        end
    end
    
    % store found Con, K values
    K_all_by_network{network_idx} = K_all_temp;
    Con_all_by_network{network_idx} = Con_all_temp;
end
%% Save analyzed data
%
%save_folder = 'H:\My Documents\Multicellular automaton\figures\two_signals\batch_sim_all_topologies_run2\count_TW';
fname_str = 'batch_sim_all_topologies_run2_count_TW_analyzed_v2_by_pset';
save( fullfile(save_folder, fname_str), 'networks_sel', 'TW_count_strict',...
    'TW_count_loose', 'K_all_by_network', 'Con_all_by_network');
%}
%% Load analyzed data
%save_folder = 'H:\My Documents\Multicellular automaton\figures\two_signals\batch_sim_all_topologies_run2\count_TW';
fname_str = 'batch_sim_all_topologies_run2_count_TW_analyzed_v2_by_pset';
load( fullfile(save_folder, fname_str), 'networks_sel', 'TW_count_strict',...
    'TW_count_loose', 'K_all_by_network', 'Con_all_by_network');

%% Plot unnormalized results
networks_idx2 = 1:5; %[1 3 6 7 8]; %[15 19 33 34 36];
num_params = [5 5 6 6 6]; % number of parameters per network

% get # parameter sets that gave at least one TW
y_data = zeros(5, 1);
for i=1:5
    idx_temp = find( TW_count_loose(i,:)==1 );
    y_data(i) = numel(idx_temp);
end
y_data = y_data/n_pset;

% overall fraction of simulations that gave TWs
%y_data = TW_count_loose(networks_idx2)/(n_pset*nsim);

h = figure;
bar( 1:numel(networks_idx2), y_data);
%set(gca, 'XTickLabel', sprintfc('Network %d', networks_sel(networks_idx2)) );
xlabel('Network');
set(gca, 'XTickLabel', sprintfc('%d', networks_sel(networks_idx2)) );
xtickangle(45)
%ylabel('Count');
ylabel('Frequency');
%set(gca, 'YScale', 'log');
set(gca, 'FontSize', 32);
%ylim([10^(-4) 1]);
ylim([0 2*10^(-3)]);
box on

% Save figures
set(h, 'Units', 'Inches', 'Position', [0.1 0.1 10 8]);
qsave = 0;
fname_str = 'Robustness_frac_sim_TW_raw_v2';
save_figure(h, 10, 8, fullfile(save_folder, fname_str),'.pdf', qsave);

%% Plot normalized results
% data
x_data = 1:numel(networks_idx2);
%y_data = (TW_count_loose(networks_idx2)'/(n_pset*nsim)).^(1./num_params);
y_data_norm = y_data.^(1./num_params');

% plot
h = figure;
bar( x_data, y_data_norm );
set(gca, 'XTickLabel', sprintfc('%d', networks_sel(networks_idx2)) );
set(gca, 'FontSize', 32);
xtickangle(45)
%ylabel('Count');
xlabel('Network');
ylabel('Normalized frequency');
%set(gca, 'YScale', 'log');
%ylim([10^(-4) 1]);
ylim([0 1]);
text(x_data-0.4, y_data_norm+0.05, sprintfc('$n_P = %d$', num_params),...
    'FontSize', 24, 'interpreter', 'latex' );

% Save figures
set(h, 'Units', 'Inches', 'Position', [0.1 0.1 10 8]);
qsave = 0;
fname_str = 'Robustness_frac_sim_TW_normalized_v2';
save_figure(h, 10, 8, fullfile(save_folder, fname_str),'.pdf', qsave);

%% Plot parameter sets as spider plots
for idx_loop=1:numel(networks_sel)
    network = networks_sel(idx_loop);
    
    K_all_temp = K_all_by_network{idx_loop};
    Con_all_temp = Con_all_by_network{idx_loop};
    
    % filter data based on present interactions (not generalized)
    if idx_loop<3
        K_idx = 1:3;
    else
        K_idx = 1:4;
    end
    
    % input vars
    %P_data = log10([Con_all_temp, K_all_temp(:,K_idx)]);
    P_data = [Con_all_temp, K_all_temp(:,K_idx)];
    %P_labels = {'$C_{ON}^{(1)}$', '$C_{ON}^{(2)}$', '$K^{(11)}$',...
    %	'$K^{(12)}$', '$K^{(21)}$', '$K^{(22)}$'};
    P_labels = {'C_{ON}^{(1)}', 'C_{ON}^{(2)}', 'K^{(11)}',...
        'K^{(12)}', 'K^{(21)}', 'K^{(22)}'};
    axes_interval = 3;
    
    % plot
    %{
    spider_plot(P_data, P_labels([1:2 K_idx+2]), axes_interval,...
        'Marker', 'o',...
        'LineStyle', '-',...
        'Color', [1 0 0],...
        'LineWidth', 2,...
        'MarkerSize', 2);
    %}
    spider_plot_linear(P_data, P_labels([1:2 K_idx+2]), axes_interval);
    %title(sprintf('n=%d, nw %d', size(P_data,1), network),...
    %    'Fontweight', 'bold', 'FontSize', 28);
    set(gcf, 'Units', 'Inches', 'Position', [2 2 10 8]);
    h = gcf;
    hold off 
    
    % save plot
    set(h, 'Units', 'Inches', 'Position', [0.1 0.1 10 8]);
    qsave = 1;
    fname_str = sprintf('Parameters_TW_spider_plot_network_%d_n_%d_linear_filled_v2',...
        network, size(P_data,1));
    save_figure(h, 10, 8, fullfile(save_folder, fname_str),'.pdf', qsave);
    pause(0.2);
    close all
end