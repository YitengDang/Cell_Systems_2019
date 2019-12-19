%% Plots trajectories with almost similar starting conditions and their Hamming distance
% (1) Hamiltonian map with p,I trajectories
% (2) Hamiltonian trajectories
% (3) Hamming distance trajectories
% (4) Final Hamming distance distribution
% lattice parameters
clear all
close all
clc
%%
gridsize = 11;
N = gridsize^2;
a0 = 5.4;
Rcell = 0.2*a0;
K = 8;
Con = 16;
p0 = 0.3;
iniON = round(p0*N);
hill = 2;
prec = 8;
%ini_noise = 0.1;
ini_noise_all = 10.^(-2:0.2:0);
n_trials = 100;

% Filename for saving
fname_out = strrep(sprintf('N%d_a0_%.1f_K_%d_Con_%d_n%d_hill%.2f', ...
    N, a0, K, Con, iniON, hill), '.', 'p'); 

% Path to search for the saved data. It searchs by the name, defined by the
% parameters chosen 
load_folder = 'H:\My Documents\Multicellular automaton\figures\finite_Hill\perturb_ini\data_binary_constraint_v2\vs_perturbation_strength';

% Calculate fN
[dist, pos] = init_dist_hex(gridsize, gridsize);
dist_vec = dist(1,:);
r = a0*dist_vec(dist_vec>0); % exclude self influence
fN = sum(Rcell*sum(exp(Rcell-r)./r)); % calculate signaling strength

% Find uniform lattice fixed points
fp = zeros(3, 1);
x0 = [0.03 0.2 0.65]; %estimates based on previous graph
%hfunc = @update_function_uniform;
hfunc2 = @(x) update_function_uniform(x, hill, Con, K, fN) - x; % update function(x) = x
for idx=1:3
    fp(idx) = fzero(hfunc2, x0(idx));
end

% Get all file names in the directory
listing = dir(load_folder);
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
%% Calculate all distances
distances_all = zeros( numel(ini_noise_all), n_trials );
for loop_idx = 1:numel(ini_noise_all)
    ini_noise = ini_noise_all(loop_idx);
    %ini_noise = 0.1;
    
    straux = '(\d+)';
    fpattern = strrep(sprintf('perturb_ini_N%d_pini%.3f_a0_%.1f_K_%d_Con_%d_hill_%.2f_ini_noise_%.2f-v%s', ...
	N, p0, a0, K, Con, hill, ini_noise, straux), '.', 'p');
    %% Load base trajectory (used for flipping)
    for i = 1:numel(names)
        fpattern0 = strrep(sprintf('perturb_ini_N%d_pini%.3f_a0_%.1f_K_%d_Con_%d_hill_%.2f_ini_noise_%.2f-v0', ...
                N, p0, a0, K, Con, hill, ini_noise), '.', 'p');
        [tokens, ~] = regexp(names{i},fpattern0,'tokens','match');
        if numel(tokens)>0
            disp('yes');
            load(fullfile(load_folder,strcat(names{i},'.mat')), 'cells_hist', 'fN');
            cells_base = cells_hist;
            cells_base_out_binary = cells_hist{end} > fp(2);
        end
    end
    %% Load other trajectories and check Hamming distance
    count = 0;
    for i = 1:numel(names)
        [tokens, ~] = regexp(names{i}, fpattern, 'tokens', 'match');
        if numel(tokens)>0
            disp(count+1);
            v = str2double(tokens{1}{1}); 
            if v>0
                count = count + 1;
                load(fullfile(load_folder, strcat(names{i},'.mat')), 'cells_hist', 'fN');
                
                % calculate distance from unflipped trajectory
                cells_perturbed_out_binary = cells_hist{end} > fp(2);
                distances_all(loop_idx, count) = sum(abs(cells_base_out_binary-cells_perturbed_out_binary));
            end
        end
    end 
end

%% Save analyzed data

fname_str = 'analyzed_vs_perturbation_strength';
fname = fullfile(load_folder, fname_str);
save(fname, 'distances_all', 'N', 'a0', 'K', 'Con', 'hill', 'p0', 'Rcell');

%% Plot Hamming distance: mean
x_data = ini_noise_all;
y_data = mean( distances_all, 2 )/N;
sigma_data = std( distances_all, 1, 2)/N;

h = figure;
errorbar(x_data, y_data, sigma_data, 'LineWidth', 2);
set(gca, 'XScale', 'log');
ylim([0 1]);
xlabel('\sigma_{pert}');
ylabel('d_H/N');
set(gca, 'FontSize', 32);

qsave = 1;
save_fig_folder = 'H:\My Documents\Multicellular automaton\figures\finite_Hill\perturb_ini\binary_v2_example';
fname_str = strrep(sprintf('N%d_a0_%.1f_K_%d_Con_%d_n%d_hill%.2f_ntrials_%d_Hamming_mean_std_v2', ...
    N, a0, K, Con, iniON, hill, n_trials), '.', 'p'); 
fname = fullfile(save_fig_folder, fname_str);
save_figure(h, 10, 8, fname, '.pdf', qsave);
%% Plot Hamming distance: distribution
bins = 0:1/N:1;
bincenters = (bins(1:end-1)+bins(2:end))/2;
h_counts = zeros( numel(ini_noise_all), numel(bins)-1 );
for i=1:numel(ini_noise_all)
    h_counts(i,:) = histcounts(distances_all(i, :)/N, bins);
end
%%
h = figure;
h_im = imagesc(log10(ini_noise_all), bincenters, h_counts'/n_trials);
set(gca, 'YDir', 'normal');
%set(h_im, 'AlphaData', h_counts' > 0)
%colormap('jet')
xlabel('log_{10}(\sigma_{pert})');
ylabel('d_H/N');
set(gca, 'FontSize', 32);
set(gcf, 'Units', 'Inches', 'Position', [1 1 10 8]);
ylim([0 1]);
c=colorbar;
ylabel(c, 'Probability');
caxis([0 0.65]);

qsave = 1;
save_fig_folder = 'H:\My Documents\Multicellular automaton\figures\finite_Hill\perturb_ini\binary_v2_example';
fname_str = strrep(sprintf('N%d_a0_%.1f_K_%d_Con_%d_n%d_hill%.2f_ntrials_%d_Hamming_distribution_', ...
    N, a0, K, Con, iniON, hill, n_trials), '.', 'p'); 
fname = fullfile(save_fig_folder, fname_str);
save_figure(h, 10, 8, fname, '.pdf', qsave);