%% Extends saved simulations by a longer time
close all
clear all

%% Simulation parameters
% variable to loop over
sigma_D_all = [0.02 0.05];

% number of simulations to do 
sim_count = 20;

% other settings
% InitiateI = 0; % generate lattice with input I?
tmax_new = 10^4; % max. number of time steps 

% folder to save simulations in
%parent_folder = 'W:\staff-bulk\tnw\BN\HY\Shared\Yiteng\two_signals\moving_cells'; 
parent_folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\moving_cells';
subfolder = 'one_signal_temp';
save_folder = fullfile(parent_folder, subfolder);
            
% default file name
sim_ID = 'one_signal';
%%
% Load trajectory
t_max_old = 1000;
%folder = 'W:\staff-bulk\tnw\BN\HY\Shared\Yiteng\two_signals\moving_cells\one_signal';
folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\moving_cells\one_signal';
%ii = 1;
%jj = 1;
nruns = 20;
for ii=1:numel(sigma_D_all)
    for jj=1:nruns
        sigma_D = sigma_D_all(ii);
        fname_str = strrep(sprintf('one_signal_sigma_D_%.3f_t_out_%d-v%d', ...
            sigma_D, t_max_old, jj), '.', 'p');
        load(fullfile(folder, fname_str));

        s = save_consts_struct;
        N = s.N;
        a0 = s.a0;
        K = s.K;
        Con = s.Con;
        M_int = s.M_int;
        hill = s.hill;
        noise = s.noise;
        rcell = s.rcell;
        sigma_D = s.sigma_D;
        growth_rate = s.growth_rate;
        R_division = s.R_division;

        Rcell = rcell*a0;
        prec = 8;
        gz = sqrt(N);

        %%
        display_fig = 0;
        cells = cells_hist{end};
        positions_all = positions_all(1:end-1); % remove/do not take into account last entry 
        positions = positions_all{end};

        %-------------dynamics-------------------------------------------------
        t = numel(cells_hist)-1;
        period = Inf; %default values
        t_onset = Inf; 

        % create figure
        if display_fig
            hin = figure;
            %set(hin, 'Position', [100 100 800 800]);
            [h_cells, h_borders] = reset_cell_figure(hin, positions, rcell);
            disp_mol = 1;
            showI = 0; 
            disp(size(cells));
            update_figure_periodic_scatter(h_cells, cells, t, disp_mol, showI, a0, distances);
        end

        %[cellsOut, ~] = update_cells_two_signals_multiply_finite_Hill(cells, distances, M_int, a0,...
        % 	Rcell, Con, Coff, K, lambda, hill, noise);
        [cellsOut, ~] = ...
            update_cells_noise_hill(cells, distances, M_int, Con, K, a0, Rcell, noise, hill, prec);

        % update positions
        [positions, distances, ~] = update_cell_positions(gz, rcell, positions, distances, sigma_D);

        %%
        while t<tmax_new % changed && period==Inf 
            pause(0.01);
            t = t+1;
            disp(t);
            cells = cellsOut;
            cells_hist{end+1} = cells; %{cells(:, 1), cells(:, 2)};
            positions_all{end+1} = positions;

            % Update figure
            if display_fig
                update_figure_periodic_cell_motion(h_cells, h_borders, cells, t, disp_mol, showI, a0, distances, positions);
            end

            %[cellsOut, ~] = update_cells_two_signals_multiply_finite_Hill(cells, distances, M_int, a0,...
            %    Rcell, Con, Coff, K, lambda, hill, noise);
            [cellsOut, ~] = ...
                update_cells_noise_hill(cells, distances, M_int, Con, K, a0, Rcell, noise, hill, prec);

            % update positions
            [positions, distances, ~] = update_cell_positions(gz, rcell, positions, distances, sigma_D);
        end
        t_out = t; %default t_out

        fprintf('Final: t_out = %d, period %d \n', t_out, period);
        %----------------------------------------------------------------------
        %% Save result
        % Format of output filename
        fname_str = strrep(sprintf('%s_sigma_D_%.3f_t_out_%d-v%d',...
            sim_ID, sigma_D, t_out, jj), '.','p');
        ext = '.mat';
        fname = fullfile(save_folder, strcat(fname_str, ext) );
        % check if filename already exists

        %{
        i=1;
        fname = fullfile(save_folder, strcat(fname_str, '-v', num2str(i), ext));
        while exist(fname, 'file') == 2
            i=i+1;
            fname = fullfile(save_folder, strcat(fname_str, '-v', num2str(i), ext));
        end
        %}

        % save
        save(fname, 'save_consts_struct', 'cells_hist', 'positions_all', 't_out',...
            'positions', 'distances');
        fprintf('Saved simulation: %s ; \n', fname);
    end
end