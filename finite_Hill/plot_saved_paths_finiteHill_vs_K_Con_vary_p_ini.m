%% Study the spatial order of patterns in the finite Hill model
% Plot how mean final p, I change with K, Con while the rest is fixed
close all
clear all
warning off
set(0, 'defaulttextinterpreter', 'tex');
%%
% Lattice parameters
gz = 15;
N = gz^2;
a0 = 6;
rcell = 0.2;
Rcell = rcell*a0;

% loop parameters
Con_list = 1:30; %5:5:30;
K_list = 1:30; %5:5:40; 
p_ini_list = 0:0.1:1;

% circuit parameters
hill = 2; % Hill coefficient
noise = 0;

% simulation parameters
num_sims = 10;
initialID = 'uniform';
sim_ID = 'one_signal_finite_hill';
tmax = 10^3;
prec = 8;
InitiateI = 0;
p_ini = Inf;
I_ini = Inf;
cells_ini = [];

% save folder (parent)
parent_folder = fullfile('N:\tnw\BN\HY\Shared\Yiteng\one_signal_finite_Hill\dynamics',...
    '2019-04-26_vs_K_Con_a0_0p50_1p50_hill2p00_vary_p_ini'); 

% Save figures folder
save_fig_folder = 'H:\My Documents\Thesis\Single gene finite Hill\Fig3';

% Calculate fN (check)
%{
mcsteps = 0;
nodisplay = 1;
[pos, dist] = initial_cells_random_markov_periodic(gz, mcsteps, rcell, nodisplay);
dist_vec = a0*dist(1,:);
r = dist_vec(dist_vec>0); % exclude self influence
fN = sum(sinh(Rcell)*sum(exp(Rcell-r)./r)); % calculate signaling strength
%}
%% Load data
% variables to store
p_final_all = zeros(numel(K_list), numel(Con_list), numel(p_ini_list), num_sims);
I_final_all = p_final_all;
t_out_all = p_final_all;
sigma_all = p_final_all;

for idx_K = 1:numel(K_list)
    this_K = K_list(idx_K);
    for idx_Con = 1:numel(Con_list)
        this_Con = Con_list(idx_Con);
        for idx_p_ini = 1:numel(p_ini_list)
            p_ini = p_ini_list(idx_p_ini);
            
            % filename pattern
            fpattern = strrep(sprintf(...
                'N%d_a0_%.2f_K_%d_Con_%d_hill_%.2f_p_ini_%.2f_t_out_%s-v%s',...
                N, a0, this_K, this_Con, hill, p_ini,...
                '(\d+)', '(\d+)'), '.', 'p');
            
            % Get all file names in the directory
            subfolder = strrep(sprintf('a0_%.2f', a0), '.', 'p');
            subfolder2 = sprintf('Con_%d_K_%d', this_Con, this_K);
            load_folder = fullfile(parent_folder, subfolder, subfolder2);
            
            listing = dir(load_folder);
            num_files = numel(listing)-2; %first two entries are not useful
            names = {};
            count = 0;
            for i = 1:num_files
                filename = listing(i+2).name;
                % remove extension and do not include txt files
                [~,name,ext] = fileparts(filename);
                if strcmp(ext, '.mat') && ~isempty(regexp(name, fpattern, 'match'))
                    count = count + 1;
                    names{count} = name;
                end
            end
            
            if numel(names)~=num_sims
                error('Wrong number of simulations in folder!');
            end
            
            % load files
            for i = 1:numel(names)
                % first get the filename given by the student
                [tokens, ~] = regexp(names{i}, fpattern, 'tokens', 'match');
                if numel(tokens) > 0
                    disp(names{i}) % displays the file name

                    % load the data
                    load(fullfile(load_folder, strcat(names{i},'.mat')),...
                        'cells_hist', 'distances', 't_out');

                    % store relevant data
                    p_final_all(idx_K, idx_Con, idx_p_ini, i) = sum(cells_hist{end})/N;
                    I_final_all(idx_K, idx_Con, idx_p_ini, i) = moranI(cells_hist{end}, a0*distances);
                    t_out_all(idx_K, idx_Con, idx_p_ini, i) = t_out;
                    sigma_all(idx_K, idx_Con, idx_p_ini, i) = std(cells_hist{end});
                end
            end
        end
    end
end

% Save loaded data
fname_str = strrep(sprintf('analyzed_data_gz_%d_a0_%.2f_K_%d_to_%d_Con_%d_to_%d_p_ini_%.2f_to_%.2f',...
    gz, a0, K_list(1), K_list(end), Con_list(1), Con_list(end),...
    p_ini_list(1), p_ini_list(end)), '.', 'p');
fname = fullfile(parent_folder, fname_str);
save( fname, 'p_final_all', 'I_final_all', 't_out_all', 'sigma_all');
%}
%% Load saved data
fname_str = strrep(sprintf('analyzed_data_gz_%d_a0_%.2f_K_%d_to_%d_Con_%d_to_%d_p_ini_%.2f_to_%.2f',...
    gz, a0, K_list(1), K_list(end), Con_list(1), Con_list(end),...
    p_ini_list(1), p_ini_list(end)), '.', 'p');
fname = fullfile(parent_folder, fname_str);
load( fname, 'p_final_all', 'I_final_all', 't_out_all', 'sigma_all');

%% Process data for scatter plot
x_data = zeros(size(p_final_all));
%x_data = zeros(30, 30, 9, 10);
y_data = x_data;
for i=1:numel(K_list)
    x_data(i,:,:,:) = K_list(i);
end
for j=1:numel(Con_list)
    y_data(:,j,:,:) = Con_list(j);
end
%p_data = p_final_all(1:30, :, :, :);
%f = fit([x_data(:), y_data(:)], p_final_all(:), 'Linear');

%% Scatter plot with all final <X_i> 
h=figure;
%plot(f);
hold on
sz = 50;
%scatter3(x_data(:), y_data(:), p_final_all(:), sz, 'filled');
%scatter3(x_data(:), y_data(:), p_final_all(:), sz, p_final_all(:), 'filled');
p=scatter3(y_data(:), x_data(:), p_final_all(:), sz, p_final_all(:), 'filled');
%scatter3(y_data(:), x_data(:), p_data(:), sz, p_data(:), 'filled');
p.MarkerFaceAlpha = 0.7;
set(gca, 'YDir', 'normal')
%xlabel('K');
%ylabel('C_{ON}');
ylabel('K');
xlabel('C_{ON}');
zlabel('Final \langle X_k \rangle');
set(gca, 'FontSize', 32);
set(h, 'Units', 'inches', 'Position', [2 2 10 8])
%c = colorbar;
%caxis([0 1]);
view([-130 30]);
pbaspect([1 1 0.5])

qsave = 0;
fname_str = strrep(sprintf(...
    'simulated_all_final_p_vary_p_gz_%d_a0_%.2f_K_%d_to_%d_Con_%d_to_%d_scatter3',...
    gz, a0, K_list(1), K_list(end), Con_list(1), Con_list(end)), '.', 'p');
path_out = fullfile(save_fig_folder, fname_str);
save_figure(h, 0, 0, path_out, '.pdf', qsave)
%{
K_idx = 12; 
Con_idx = 30;
squeeze(p_final_all(K_idx, Con_idx, :, :))
%}
%% Obtain number of SS values (only for 1 or 2)
% Classify data into binary classes (high/low final p) 
p_min_all = min(p_final_all(:,:,:), [], 3);
p_max_all = max(p_final_all(:,:,:), [], 3);
two_ss = (p_max_all-p_min_all>0.0001 ); % two distinct SS if difference between values is sufficiently large

h=figure;
imagesc(K_list, Con_list, two_ss' );
%imagesc( Con_list, K_list, two_ss );
%colorbar;
set(gca, 'YDir', 'normal');
xlabel('K');
ylabel('C_{ON}');
set(gca, 'FontSize', 32);
set(gca, 'XTick', [1 5:5:30], 'YTick', [1 5:5:30]);
set(h, 'Units', 'inches', 'Position', [2 2 10 8])
pbaspect([1 1 1]);

qsave = 0;
fname_str = strrep(sprintf(...
    'simulated_num_SS_vary_p_gz_%d_a0_%.2f_K_%d_to_%d_Con_%d_to_%d_imagesc',...
    gz, a0, K_list(1), K_list(end), Con_list(1), Con_list(end)), '.', 'p');
path_out = fullfile(save_fig_folder, fname_str);
save_figure(h, 0, 0, path_out, '.pdf', qsave)

% Check that the variability within each class is small -> DONE
%{
p_data_binary = zeros( size(p_final_all) ); % Binary data
p_data_var = zeros( numel(K_list), numel(Con_list), 2 ); % variability in p_data after binary classification
for i=1:size(p_final_all,1)
    for j=1:size(p_final_all,2)
        this_p_min = p_min_all(i,j);
        this_p_max = p_max_all(i,j);
        this_p_data = squeeze(p_final_all(i,j,:,:));
        % compare whether the data is closer to the minimum value or the
        % maximum value
        this_p_data_binary = (abs( this_p_data - this_p_min ) > abs( this_p_data - this_p_max ));
        p_data_binary(i,j,:,:) = this_p_data_binary;
        p_data_var(i, j, 1) = std( this_p_data(this_p_data_binary) );
        p_data_var(i, j, 2) = std( this_p_data(~this_p_data_binary) );
    end
end

figure;
subplot(2,1,1);
imagesc(p_data_var(:,:,1));
colorbar;
subplot(2,1,2);
imagesc(p_data_var(:,:,2));
colorbar;
%}
%% Plot SEM of <X_i>(t_final) vs K and Con
p_end_SEM = std(mean(p_final_all, 4), 1, 3)/sqrt(numel(p_ini_list));

h=figure;
imagesc(K_list, Con_list, p_end_SEM')
set(gca, 'YDir', 'normal')
ylabel('C_{ON}');
xlabel('K');
c=colorbar;
%caxis([0 1]);
set(gca, 'FontSize', 32);
ylabel(c, 'Final gene expression SEM');
set(h, 'Units', 'inches', 'Position', [2 2 10 8])
pbaspect([1 1 1]);

qsave = 0;
fname_str = strrep(sprintf('simulated_mean_final_p_SEM_vary_p_gz_%d_a0_%.2f_K_%d_to_%d_Con_%d_to_%d_heatmap',...
    gz, a0, K_list(1), K_list(end), Con_list(1), Con_list(end)), '.', 'p');
path_out = fullfile(save_fig_folder, fname_str);
save_figure(h, 0, 0, path_out, '.pdf', qsave)

%% Plot sigma(X_i(t_final)) vs K and Con
sigma_mean = mean(sigma_all(:,:,:), 3);

h=figure;
imagesc(K_list, Con_list, sigma_mean')
set(gca, 'YDir', 'normal')
ylabel('C_{ON}');
xlabel('K');
c=colorbar;
%caxis([0 0.3]);
set(gca, 'FontSize', 32);
ylabel(c, 'Mean variability of final state');
set(h, 'Units', 'inches', 'Position', [2 2 10 8])
pbaspect([1 1 1]);

qsave = 0;
fname_str = strrep(sprintf('simulated_mean_sigma_final_vary_p_gz_%d_a0_%.2f_K_%d_to_%d_Con_%d_to_%d_heatmap',...
    gz, a0, K_list(1), K_list(end), Con_list(1), Con_list(end)), '.', 'p');
path_out = fullfile(save_fig_folder, fname_str);
save_figure(h, 0, 0, path_out, '.pdf', qsave)
%% Plot homogeneous vs heterogeneous final state vs K and Con
% Is it possible to have a heterogeneous final state?
sigma_max = max(sigma_all(:,:,:),[], 3);
sigma_threshold = 0.01;

h=figure;
imagesc(K_list, Con_list, sigma_max' > sigma_threshold)
set(gca, 'YDir', 'normal')
ylabel('C_{ON}');
xlabel('K');
colormap('winter')
%c=colorbar;
%caxis([0 0.3]);
set(gca, 'XTick', [1 5:5:30]);
set(gca, 'FontSize', 32);
%ylabel(c, 'Mean variability of final state');
set(h, 'Units', 'inches', 'Position', [2 2 10 8])
pbaspect([1 1 1]);

qsave = 0;
fname_str = strrep(sprintf('simulated_max_sigma_final_geq_%.2f_vary_p_gz_%d_a0_%.2f_K_%d_to_%d_Con_%d_to_%d_heatmap',...
    sigma_threshold, gz, a0, K_list(1), K_list(end), Con_list(1), Con_list(end)), '.', 'p');
path_out = fullfile(save_fig_folder, fname_str);
save_figure(h, 0, 0, path_out, '.pdf', qsave)

%% Plot t_out vs K and Con
t_out_mean = mean(mean(t_out_all, 3), 4);

h=figure;
imagesc(Con_list, K_list, t_out_mean)
set(gca, 'YDir', 'normal')
xlabel('C_{ON}');
ylabel('K');
c=colorbar;
caxis([0 tmax]);
set(gca, 'FontSize', 32);
ylabel(c, 'Mean simulation time');
set(h, 'Units', 'inches', 'Position', [2 2 10 8])

qsave = 0;
fname_str = strrep(sprintf('simulated_mean_t_eq_final_vary_p_gz_%d_a0_%.2f_K_%d_to_%d_Con_%d_to_%d_heatmap',...
    gz, a0, K_list(1), K_list(end), Con_list(1), Con_list(end)), '.', 'p');
path_out = fullfile(save_fig_folder, fname_str);
save_figure(h, 0, 0, path_out, '.pdf', qsave)

%% Plot <I>(t_final) vs K and Con
% N.B. I is not a reliable measure for spatial order anymore 
I_end_mean = mean(mean(I_final_all, 3), 4);

h=figure;
imagesc(Con_list, K_list, I_end_mean)
set(gca, 'YDir', 'normal')
xlabel('C_{ON}');
ylabel('K');
c=colorbar;
%caxis([0.25 0.4]);
set(gca, 'FontSize', 32);
ylabel(c, 'Mean final spatial order');
set(h, 'Units', 'inches', 'Position', [2 2 10 8])

qsave = 0;
fname_str = strrep(sprintf('simulated_mean_final_I_vary_p_gz_%d_a0_%.2f_K_%d_to_%d_Con_%d_to_%d_heatmap',...
    gz, a0, K_list(1), K_list(end), Con_list(1), Con_list(end)), '.', 'p');
path_out = fullfile(save_fig_folder, fname_str);
save_figure(h, 0, 0, path_out, '.pdf', qsave)