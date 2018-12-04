function [cells_hist, period, t_onset] = time_evolution_save_func_one_signal_moving_cells(...
    N, a0, Rcell, hill, noise, M_int, K, Con,...
    distances, positions, sigma_D, cells_ini, growth_rate, R_division, ...
    sim_ID, tmax, save_folder, display_fig)
    %%
    % Initialize
    gz = sqrt(N);
    rcell = Rcell/a0;
    cells_hist = {};
    prec = 8; % precision for finite Hill
    
    % generate initial lattice
    cells_hist{1} = cells_ini;
    positions_all = {};
    positions_all{1} = positions;
    %-------------dynamics-------------------------------------------------
    t = 0;
    period = Inf; %default values
    t_onset = Inf; 
    cells = cells_ini;
    %[cellsOut, ~] = update_cells_two_signals_multiply_finite_Hill(cells, distances, M_int, a0,...
    % 	Rcell, Con, Coff, K, lambda, hill, noise);
    [cellsOut, ~] = ...
        update_cells_noise_hill(cells, distances, M_int, Con, K, a0, Rcell, noise, hill, prec);
    
    % update positions
    [positions, distances, ~] = update_cell_positions(gz, rcell, positions, distances, sigma_D);
    
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
    % t_check = 10^3; 
    while t<tmax % changed && period==Inf 
        pause(0.01);
        t = t+1; 
        disp(t);
        cells = cellsOut;
        cells_hist{end+1} = cells; %{cells(:, 1), cells(:, 2)};
        positions_all{end+1} = positions;
        %{
        if mod(t, t_check)==0
            [period, t_onset] = periodicity_test_short(cells_hist); 
        end
        %}
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
    
    % if periodicity found, refine check to find period
    %{
    if period<Inf && t>t_ac
        [period, t_onset] = periodicity_test_detailed(cells_hist, t_check, period);
        t_out = t_onset + period;
    end
    %}
    
    fprintf('Final: t_out = %d, period %d \n', t_out, period);
    %----------------------------------------------------------------------
    % Save result       
    % Format of output filename
    fname_str = strrep(sprintf('%s_sigma_D_%.3f_t_out_%d',...
    	sim_ID, sigma_D, t_out), '.','p');
    ext = '.mat';
    
    % check if filename already exists
    i=1;
    fname = fullfile(save_folder, strcat(fname_str, '-v', num2str(i), ext));
    while exist(fname, 'file') == 2
        i=i+1;
        fname = fullfile(save_folder, strcat(fname_str, '-v', num2str(i), ext));
    end
    
    % variables to save
    rcell = Rcell/a0;
    save_vars = {N, a0, K, Con, M_int, hill,...
        noise, cells_ini, rcell, Rcell,...
        sim_ID, tmax, sigma_D, growth_rate, R_division};
    save_vars_lbl = {'N', 'a0', 'K', 'Con', 'M_int', 'hill',...
        'noise', 'cells_ini', 'rcell', 'Rcell',...
        'sim_ID', 'tmax', 'sigma_D', 'growth_rate', 'R_division'};
    save_consts_struct = cell2struct(save_vars, save_vars_lbl, 2);
    
    % save
    save(fname, 'save_consts_struct', 'cells_hist', 'positions_all', 't_out',...
        'positions', 'distances');
    fprintf('Saved simulation: %s ; \n', fname);
end