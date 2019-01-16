function [cells_hist, period, t_onset] = time_evolution_save_func_efficient_checks_moving_cells(...
    N, a0, Rcell, lambda, hill, noise, M_int, K, Con, Coff,...
    distances, positions, mcsteps, sigma_D, cells_ini, growth_rate, R_division, ...
    sim_ID, tmax, save_folder, fname_str_template, display_fig)
    % Similar to time_evolution_save_func, but
    % implements an efficient periodicity test by checking only every
    % t_check time steps whether trajectory has become periodic (after a
    % transient phase of t_ac steps when it checks every step).
    %%
    % Initialize
    gz = sqrt(N);
    rcell = Rcell/a0;
    cells_hist = {};
    
    % randomize initial lattice
    if isempty(positions) && isempty(distances)
        nodisplay = 0;
        [positions, distances] = initial_cells_random_markov_periodic(gz, mcsteps,...
            rcell, nodisplay);
    end
    
    % note: initial state given as function input, generate random state
    % outside of function
    
    % store initial state and lattice
    cells_hist{1} = cells_ini;
    positions_all = {};
    positions_all{1} = positions;
        
    
    % Show initial condiguration
    if display_fig
        hin = figure;
        %set(hin, 'Position', [100 100 800 800]);
        [h_cells, h_borders] = reset_cell_figure(hin, positions, rcell);
        disp_mol = 12;
        showI = 0; 
        t=0;
        %disp(size(cells));
        update_figure_periodic_scatter(h_cells, cells_ini, t, disp_mol,...
            showI, a0, distances);
        pause(1);
    end
    %-------------dynamics-----------------------------------------
    t = 0;
    period = Inf; %default values
    t_onset = Inf; 
    cells = cells_ini;
    
    % update cells
    [cellsOut, ~] = update_cells_two_signals_multiply_finite_Hill(cells, distances, M_int, a0,...
    	Rcell, Con, Coff, K, lambda, hill, noise);
    
    % update positions
    [positions, distances, ~] = update_cell_positions(gz, rcell, positions, distances, sigma_D);
    
    %%
    %{
    t_ac = 10^2;
    while t<t_ac % changed && period==Inf
        pause(0.01);
        t = t+1;
        cells = cellsOut;
        cells_hist{end+1} = cells; %{cells(:, 1), cells(:, 2)};
        % [period, t_onset] = periodicity_test_short(cells_hist); 
        
        % Update figure
        if display_fig
        update_figure_periodic_cell_motion(h_cells, h_borders, cells, t, disp_mol, showI, a0, distances, positions);
        end
    
        [cellsOut, ~] = update_cells_two_signals_multiply_finite_Hill(cells, distances, M_int, a0,...
            Rcell, Con, Coff, K, lambda, hill, noise);
        
        % update positions
        [positions, distances, ~] = update_cell_positions(gz, rcell, positions, distances, sigma_D);
        pos_hist{end+1} = positions;
    end
    %}
    % check periodically after t_ac time steps, with period t_check
    %t_check = 10^3; 
    
    %periodicity_over_time = []; % store all found periodicities over time
    periodicity_vs_t = []; % store all found periodicities over time
    t_onset_over_time = []; % store onset times of periodic trajectories
    %t_check = 1; % start checking periodicity from time 1
    trav_wave_2_vs_t = []; % whether simulation was a TW (according to second test, only constant p) for one period  
    
    while t<tmax % && changed %&& period==Inf
        pause(0.001);
        t = t+1; 
        %disp(t);
        cells = cellsOut;
        cells_hist{end+1} = cells; %{cells(:, 1), cells(:, 2)};
        positions_all{end+1} = positions;
        %
        
        % Check for periodicity
        [period, t_onset] = periodicity_test_short(cells_hist);
        if period<Inf
            %[period, t_onset] = periodicity_test_detailed(cells_hist, t_check, period);
            [period, t_onset] = periodicity_test_short_reversed(cells_hist(1:t+1));
            periodicity_vs_t(end+1) = period;
        else
            periodicity_vs_t(end+1) = 0;
        end
        %periodicity_over_time(end+1) = period;
        t_onset_over_time(end+1) = t_onset;
        
        % Check for travelling waves
        if period<Inf
            cells_hist_temp = cells_hist(t-period+1:t+1);
            %dist = ones(N); % not needed, random value
            [~, trav_wave_2] = travelling_wave_test(cells_hist_temp, a0,...
                period, numel(cells_hist_temp)-1 , distances);
            trav_wave_2_vs_t(end+1) = trav_wave_2;
            fprintf('Time = %d, TW? %d \n', t, trav_wave_2);
        else
            trav_wave_2_vs_t(end+1) = 0;
        end
        %}
        
        % Update figure
        if display_fig
            update_figure_periodic_cell_motion(h_cells, h_borders, cells, t, disp_mol, showI, a0, distances, positions);
        end
        
        % update cell states
        [cellsOut, ~] = update_cells_two_signals_multiply_finite_Hill(cells,...
            distances, M_int, a0,...
            Rcell, Con, Coff, K, lambda, hill, noise);
        
        % update positions
        [positions, distances, ~] = update_cell_positions(gz, rcell,...
            positions, distances, sigma_D);        
    end

    t_out = t; %default t_out
    
    % if periodicity found, refine check to find period
    %{
    if period < Inf && t > t_ac
        [period, t_onset] = periodicity_test_detailed(cells_hist, t_check, period);
        t_out = t_onset + period;
    end
    %}
    
    %fprintf('Final: t_out = %d, period %d \n', t_out, period);
    fprintf('Final: t_out = %d, periodicity found? %d \n', t_out,...
        min(periodicity_vs_t) );
    
    %----------------------------------------------------------------------
    % Save result       
    % Format of output filename
    %fname_str = strrep(sprintf('%s_sigma_D_%.3f_t_out_%d',...
    % 	sim_ID, sigma_D, t_out), '.','p');
    fname_str =  strrep(sprintf('%s_tmax_%d', fname_str_template,...
        tmax), '.','p');
    
    ext = '.mat';
    
    % check if filename already exists
    v=1; 
    fname = fullfile(save_folder, strcat(fname_str, '-v', num2str(v), ext));
    while exist(fname, 'file') == 2
        v=v+1;
        fname = fullfile(save_folder, strcat(fname_str, '-v', num2str(v), ext));
    end
    
    % variables to save
    rcell = Rcell/a0;
    lambda12 = lambda(2);
    save_vars = {N, a0, K, Con, Coff, M_int, hill,...
        noise, rcell, Rcell, lambda12, lambda,...
        sim_ID, tmax, sigma_D, growth_rate, R_division};
    save_vars_lbl = {'N', 'a0', 'K', 'Con', 'Coff', 'M_int', 'hill',...
        'noise', 'rcell', 'Rcell',  'lambda12', 'lambda', ...
        'sim_ID', 'tmax', 'sigma_D', 'growth_rate', 'R_division'};
    save_consts_struct = cell2struct(save_vars, save_vars_lbl, 2);
    
    % save
    save(fname, 'save_consts_struct', 'cells_hist', 'positions_all', 't_out',...
        'positions', 'distances',...
        'periodicity_vs_t', 'trav_wave_2_vs_t');
    
    fprintf('Saved simulation: %s ; \n', fname);
end