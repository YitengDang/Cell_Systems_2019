% Analyze batch simulation results for moving cells 
% Plot pin-pout maps (1 signal)
close all
clear all
% maxNumCompThreads(4);
% warning off
set(0, 'defaulttextinterpreter', 'latex');

%% Set folders
parent_folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\moving_cells_one_signal';

% Save data folder
save_folder = fullfile(parent_folder, 'figures');

% Save figures folder
save_fig_folder = fullfile(parent_folder, 'figures');
%% Load simulations for moving cells
M_int = 1;
gz = 15;
N = gz^2;
rcell_sigma = 0;
sigma_D = 0;
t_max = 100;
mcsteps = 10^4;

subfolder = 'random_position_non_moving';
folder = fullfile(parent_folder, subfolder);

% Load growing/moving cells simulations
pin_pout_count = zeros(N+1);

% get all file names
files = dir(folder);
files = files(3:end);
names_all = {files.name};

for ii=1:numel(names_all)
    filename = names_all{ii};
    
    %fname_str = 'one_signal_sigma_D_0p001_t_out_1000-v1';
    
    pattern = strrep(sprintf('one_signal_M_int_%d_mcsteps_%d_sigma_D_%.3f_no_growth_tmax_%d_iniON_%s_t_out_%s-v%s',...
        M_int, mcsteps, sigma_D, t_max, '(\d+)', '(\d+)', '(\d+)'), '.' ,'p');
    [tokens, ~] = regexp(  filename, pattern, 'tokens', 'match');

    if ~isempty(tokens)
        fname = fullfile(folder, filename);
        disp(fname);
        load(fname, 'cells_hist');
        %load(fname, 'save_consts_struct', 'cells_hist', 't_out',...
        %    'changed', 'positions', 'distances', 'positions_all', 'rcell_hist');

        n_in = sum(cells_hist{1});
        n_out = sum(cells_hist{end});
        
        % update pin-pout count
        pin_pout_count(n_in+1, n_out+1) = pin_pout_count(n_in+1, n_out+1) + 1;
    end
end

%% Save analyzed data
%
save_fname_str = strrep(sprintf('pin_pout_analyzed_data_%s_t_max_%d_nruns_%d_rcell_sigma_%.1f_sigma_D_%.3f',...
    subfolder, t_max, nruns, rcell_sigma, sigma_D), '.', 'p');
save_file = fullfile(parent_folder, save_fname_str);
save(save_file, 'subfolder', 'p_final_all', 'p_final_all_nc', 'period_all',...
    'rcell_sigma', 'k_growth', 'sigma_D', 't_out');
%}

%% Plot pin-pout map
% normalize by number of p_in 
pin_pout = pin_pout_count./sum(pin_pout_count, 2);
p_all = (0:N)/N;

h = figure;
im_fig = imagesc(p_all, p_all, pin_pout');
xlim([0 1]);
ylim([0 1]);
xlabel('$p_{in}$');
ylabel('$p_{out}$');
set(gca, 'Ydir', 'normal', 'FontSize', 20);
set(im_fig, 'AlphaData', pin_pout_count' > 0);

qsave = 1;
fname_str = sprintf('one_signal_pin-pout_heat_map_%s', subfolder);
path_out = fullfile(parent_folder, fname_str);
save_figure(h, 10, 8, path_out, '.pdf', qsave)
