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
a0 = 0.5;
rcell = 0.2;
Rcell = rcell*a0;

% loop parameters
Con_list = 1:30;
K_list = 1:40;
%p_ini_list = 0.1:0.1:0.9;

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
    '2019-04-25_vs_K_Con_a0_0p50_1p50_hill2p00_uniform'); 

% Save figures folder
save_fig_folder = 'H:\My Documents\Thesis\Single gene finite Hill\Fig2';
%% Load data
% variables to store
p_final_all = zeros(numel(K_list), numel(Con_list), num_sims);
I_final_all = p_final_all;
t_out_all = p_final_all;
sigma_all = p_final_all;

for idx_K = 1:numel(K_list)
    this_K = K_list(idx_K);
    for idx_Con = 1:numel(Con_list)
        this_Con = Con_list(idx_Con);
        
        % filename pattern
        fpattern = strrep(sprintf('N%d_a0_%.2f_K_%d_Con_%d_hill_%.2f_t_out_%s-v%s',...
            N, a0, this_K, this_Con, hill,...
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
            if strcmp(ext, '.mat')
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
                p_final_all(idx_K, idx_Con, i) = sum(cells_hist{end})/N;
                I_final_all(idx_K, idx_Con, i) = moranI(cells_hist{end}, a0*distances);
                t_out_all(idx_K, idx_Con, i) = t_out;
                sigma_all(idx_K, idx_Con, i) = std(cells_hist{end});
            end
        end
        
    end
end

% Save loaded data
fname_str = strrep(sprintf('analyzed_data_gz_%d_a0_%.2f_K_%d_to_%d_Con_%d_to_%d',...
    gz, a0, K_list(1), K_list(end), Con_list(1), Con_list(end)), '.', 'p');
fname = fullfile(parent_folder, fname_str);
save( fname, 'p_final_all', 'I_final_all', 't_out_all', 'sigma_all');
%% Load saved data
fname_str = strrep(sprintf('analyzed_data_gz_%d_a0_%.2f_K_%d_to_%d_Con_%d_to_%d',...
    gz, a0, K_list(1), K_list(end), Con_list(1), Con_list(end)), '.', 'p');
fname = fullfile(parent_folder, fname_str);
load( fname, 'p_final_all', 'I_final_all', 't_out_all', 'sigma_all');

%% Scatter plot with all final <X_i> 
x_data = zeros(size(p_final_all));
for i=1:numel(K_list)
    x_data(i,:,:) = K_list(i);
end
y_data = zeros(size(x_data));
for j=1:numel(Con_list)
    y_data(:,j,:) = Con_list(j);
end

f = fit([x_data(:), y_data(:)], p_final_all(:), 'Linear');

%%
h=figure;
plot(f);
hold on
sz = 36;
%scatter3(x_data(:), y_data(:), p_final_all(:), sz, 'filled');
%scatter3(x_data(:), y_data(:), p_final_all(:), sz, p_final_all(:), 'filled');
set(gca, 'YDir', 'normal')
xlabel('K');
ylabel('C_{ON}');
zlabel('Final \langle X_k \rangle');
set(gca, 'FontSize', 32);
set(h, 'Units', 'inches', 'Position', [2 2 10 8])
%c = colorbar;
%caxis([0 1]);
view([45 30]);

qsave = 1;
fname_str = strrep(sprintf('Scatter3_all_final_p_gz_%d_a0_%.2f_K_%d_to_%d_Con_%d_to_%d_v3_fit_linear_interpolation',...
    gz, a0, K_list(1), K_list(end), Con_list(1), Con_list(end)), '.', 'p');
path_out = fullfile(save_fig_folder, fname_str);
save_figure(h, 0, 0, path_out, '.pdf', qsave)
%% Plot <X_i>(t_final) vs K and Con (avg. over trajectories)
p_end_mean = mean(p_final_all, 3);

h=figure;
imagesc(Con_list, K_list, p_end_mean)
set(gca, 'YDir', 'normal')
xlabel('C_{ON}');
ylabel('K');
c=colorbar;
caxis([0 1]);
set(gca, 'FontSize', 32);
ylabel(c, 'Mean final gene expression');
set(h, 'Units', 'inches', 'Position', [2 2 10 8])

qsave = 1;
fname_str = strrep(sprintf('Heatmap_mean_final_p_gz_%d_a0_%.2f_K_%d_to_%d_Con_%d_to_%d',...
    gz, a0, K_list(1), K_list(end), Con_list(1), Con_list(end)), '.', 'p');
path_out = fullfile(save_fig_folder, fname_str);
save_figure(h, 0, 0, path_out, '.pdf', qsave)

%% Plot SEM of <X_i>(t_final) vs K and Con
p_end_SEM = std(p_final_all, 1, 3)/sqrt(num_sims);

h=figure;
imagesc(Con_list, K_list, p_end_SEM)
set(gca, 'YDir', 'normal')
xlabel('C_{ON}');
ylabel('K');
c=colorbar;
caxis([0 1]);
set(gca, 'FontSize', 32);
ylabel(c, 'Mean final gene expression');
set(h, 'Units', 'inches', 'Position', [2 2 10 8])

qsave = 1;
fname_str = strrep(sprintf('Heatmap_mean_final_p_SEM_gz_%d_a0_%.2f_K_%d_to_%d_Con_%d_to_%d',...
    gz, a0, K_list(1), K_list(end), Con_list(1), Con_list(end)), '.', 'p');
path_out = fullfile(save_fig_folder, fname_str);
save_figure(h, 0, 0, path_out, '.pdf', qsave)

%% Plot sigma(X_i(t_final)) vs K and Con
sigma_mean = mean(sigma_all, 3);

h=figure;
imagesc(Con_list, K_list, sigma_mean)
set(gca, 'YDir', 'normal')
xlabel('C_{ON}');
ylabel('K');
c=colorbar;
caxis([0 0.3]);
set(gca, 'FontSize', 32);
ylabel(c, 'Mean variability of final state');
set(h, 'Units', 'inches', 'Position', [2 2 10 8])

qsave = 1;
fname_str = strrep(sprintf('Heatmap_mean_sigma_final_gz_%d_a0_%.2f_K_%d_to_%d_Con_%d_to_%d',...
    gz, a0, K_list(1), K_list(end), Con_list(1), Con_list(end)), '.', 'p');
path_out = fullfile(save_fig_folder, fname_str);
save_figure(h, 0, 0, path_out, '.pdf', qsave)

%% Plot t_out vs K and Con
t_out_mean = mean(t_out_all, 3);

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

qsave = 1;
fname_str = strrep(sprintf('Heatmap_mean_t_eq_final_gz_%d_a0_%.2f_K_%d_to_%d_Con_%d_to_%d',...
    gz, a0, K_list(1), K_list(end), Con_list(1), Con_list(end)), '.', 'p');
path_out = fullfile(save_fig_folder, fname_str);
save_figure(h, 0, 0, path_out, '.pdf', qsave)

%% Plot <I>(t_final) vs K and Con
% N.B. I is not a reliable measure for spatial order anymore 
I_end_mean = mean(I_final_all, 3);

h=figure;
imagesc(Con_list, K_list, I_end_mean)
set(gca, 'YDir', 'normal')
xlabel('C_{ON}');
ylabel('K');
c=colorbar;
caxis([0.25 0.4]);
set(gca, 'FontSize', 32);
ylabel(c, 'Mean final spatial order');
set(h, 'Units', 'inches', 'Position', [2 2 10 8])

qsave = 1;
fname_str = strrep(sprintf('Heatmap_mean_final_I_gz_%d_a0_%.2f_K_%d_to_%d_Con_%d_to_%d',...
    gz, a0, K_list(1), K_list(end), Con_list(1), Con_list(end)), '.', 'p');
path_out = fullfile(save_fig_folder, fname_str);
save_figure(h, 0, 0, path_out, '.pdf', qsave)