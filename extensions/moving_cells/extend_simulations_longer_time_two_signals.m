%% Extends saved simulations by a longer time, two signals
close all
clear all

%% Simulation parameters
% variable to loop over
sigma_D_all = 10.^(-4:0.5:-2);

% other variables
gz = 15;
N = gz^2;

% number of simulations to do 
sim_count = 1;

% other settings
% InitiateI = 0; % generate lattice with input I?
tmax_new = 2000; %2*10^3; % max. number of time steps 

% parameter set of loaded data
pset = 20;

% folder to save simulations in
save_folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\moving_cells_TW\TW_formation_network_15\Parameter_set_20'; 
%save_folder = fullfile(parent_folder, subfolder);
            
% default file name
sim_ID = 'two_signal_mult';

%%
% Load trajectory
t_max_old = 1000;
%folder = 'W:\staff-bulk\tnw\BN\HY\Shared\Yiteng\two_signals\moving_cells\one_signal';
%folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\moving_cells\one_signal_temp';
folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\moving_cells_TW\TW_formation_network_15\Parameter_set_20_tmax_1000';
nruns = 100;
for ii=1:numel(sigma_D_all)
    sigma_D = sigma_D_all(ii);
    %disp(sigma_D);
    for idx_sim=1:nruns
        % -------------- ad hoc fix for misnaming files ----------
        if sigma_D==10^(-3.5)
            version = 2;
        else
            version = 1;
        end
        %---------------------------------------------------------
        fname_str = strrep(sprintf('two_signal_mult_N%d_params_%d_sim_%d_sigma_D_%.3f_tmax_%d-v%d', ...
            N, pset, idx_sim, sigma_D, t_max_old, version), '.', 'p');
        disp(fname_str);
        load(fullfile(folder, fname_str), 'save_consts_struct', 'cells_hist',...
            'positions_all', 'distances');
        %
        s = save_consts_struct;
        N = s.N;
        a0 = s.a0;
        K = s.K;
        Con = s.Con;
        Coff = s.Coff;
        M_int = s.M_int;
        hill = s.hill;
        noise = s.noise;
        rcell = s.rcell;
        lambda12 = s.lambda12;
        sigma_D = s.sigma_D;
        growth_rate = s.growth_rate;
        R_division = s.R_division;
        
        Rcell = rcell*a0;
        lambda = [1 lambda12];
        prec = 8;
        gz = sqrt(N);

        %%
        display_fig = 0;
        cells = cells_hist{end};
        % only for earlier batch of simulations: correct positions_all
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
            disp_mol = 12;
            showI = 0; 
            %disp(size(cells));
            update_figure_periodic_scatter(h_cells, cells, t, disp_mol, showI, a0, distances);
        end

        [cellsOut, ~] = update_cells_two_signals_multiply_finite_Hill(cells, distances, M_int, a0,...
         	Rcell, Con, Coff, K, lambda, hill, noise);
        %[cellsOut, ~] = ...
        %    update_cells_noise_hill(cells, distances, M_int, Con, K, a0, Rcell, noise, hill, prec);

        % update positions
        [positions, distances, ~] = update_cell_positions(gz, rcell, positions, distances, sigma_D);

        %%
        while t<tmax_new % changed && period==Inf 
            t = t+1;
            disp(t);
            cells = cellsOut;
            cells_hist{end+1} = cells; %{cells(:, 1), cells(:, 2)};
            positions_all{end+1} = positions;

            % Update figure
            if display_fig
                pause(0.01);
                update_figure_periodic_cell_motion(h_cells, h_borders, cells, t, disp_mol, showI, a0, distances, positions);
            end

            [cellsOut, ~] = update_cells_two_signals_multiply_finite_Hill(cells, distances, M_int, a0,...
                Rcell, Con, Coff, K, lambda, hill, noise);
            %[cellsOut, ~] = ...
            %    update_cells_noise_hill(cells, distances, M_int, Con, K, a0, Rcell, noise, hill, prec);

            % update positions
            [positions, distances, ~] = update_cell_positions(gz, rcell, positions, distances, sigma_D);
        end
        t_out = t; %default t_out

        fprintf('Final: t_out = %d, period %d \n', t_out, period);
        %}
        %----------------------------------------------------------------------
        %% Save result
        % Format of output filename
        fname_str = strrep(sprintf('%s_N%d_sim_%d_params_%d_sigma_D_%.4f_tmax_%d',...
            sim_ID, N, pset, idx_sim, sigma_D, tmax_new), '.','p');
        ext = '.mat';
        %fname = fullfile(save_folder, strcat(fname_str, ext) );
        
        i=1;
        fname = fullfile(save_folder, strcat(fname_str, '-v', num2str(i), ext));
        % check if filename already exists
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