% Analyze batch simulation results for moving cells
close all
clear all
% maxNumCompThreads(4);
% warning off
set(0, 'defaulttextinterpreter', 'latex');

%% Set folders
parent_folder = 'H:\My Documents\Multicellular automaton\temp';

% Save data folder
save_folder = fullfile(parent_folder, 'figures');

% Save figures folder
save_fig_folder = fullfile(parent_folder, 'figures');
%% One signal
%
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
nruns = 38;
p_final_all = zeros(nruns, 2);
for jj=1:nruns
    %fname_str = 'one_signal_sigma_D_0p001_t_out_1000-v1';
    fname_str = strrep(sprintf('two_signals_rcell_sigma_%.1f_K_growth_%.1f_sigma_D_%.3f_t_out_%d-v%d',...
        rcell_sigma, k_growth, sigma_D, t_out, jj), '.' ,'p');
    fname = fullfile(folder, strcat(fname_str, '.mat'));
    if exist(fname, 'file')==2
        disp(fname);
        load(fname)
        
        cells_final = cells_hist{end};
        p_final_all(jj, :) = sum(cells_final)/N;       
    end
end

%% Load negative control simulations
files = dir(folder);
files = files(3:end);
for ii=1:numel(files)
    filename = files(ii).name;
    pattern = strrep(sprintf('two_signals_rcell_sigma_%.1f_K_growth_%.1f_sigma_D_%.3f_t_out_%s_neg_control-v%s',...
        rcell_sigma, k_growth, sigma_D, '(\d+)', '(\d+)'), '.', 'p');
    [tokens, ~] = regexp(filename, pattern, 'tokens', 'match');
    %disp(filename);
    %
    if ~isempty(tokens)
        disp(filename);
    end
    %}
end
%% Save analyzed data
%{
save_fname_str = sprintf('analyzed_data_%s_t_max_%d_nruns_%d_rcell_sigma_%.1f_K_growth_%.1f_sigma_D_%.3f',...
    subfolder, t_max, nruns, rcell_sigma, K_growth, sigma_D);
save_file = fullfile(parent_folder, save_fname_str);
save(save_file, 'subfolder', 'p_final_all', 'rcell_sigma', 'K_growth',...
    'sigma_D', 't_out');
%}

%% Plot final p
bins = 0:0.1:1;
h=figure;
hold on
histogram( p_final_all(:,1), bins );
histogram( p_final_all(:,2), bins );
title(sprintf('n=%d', nruns));
xlabel('$p_{final}$');
ylabel('Count');
set(gca, 'FontSIze', 20);