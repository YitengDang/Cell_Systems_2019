function [cells_hist, period, t_onset] = check_TW_propagation_save(...
    N, a0, Rcell, lambda, hill, noise, M_int, K, Con, Coff,...
    distances, positions, sim_ID, mcsteps, InitiateI, p0, I0, cells_ini, tmax,...
    save_folder, fname_str_template, display_fig)
    % Runs the simulation only for gz time steps and checks whether a TW is
    % obtained again
    
    cells_hist = {};
    gz = sqrt(N);
    rcell = Rcell/a0;
    
    % store initial config
    cells = cells_ini;
    cells_hist{end+1} = cells; %{cells(:, 1), cells(:, 2)};
    
    % Show initial condiguration
    if display_fig
        hin = figure;
        %set(hin, 'Position', [100 100 800 800]);
        [h_cells, h_borders] = reset_cell_figure(hin, positions, rcell);
        disp_mol = 12;
        showI = 0; 
        t=0;
        %disp(size(cells));
        update_figure_periodic_scatter(h_cells, cells, t, disp_mol, showI, a0, distances);
        pause(1);
    end
    
    %-------------dynamics-------------------------------------------------
    t = 0;
    period = Inf; %default values
    t_onset = Inf; 
    [cellsOut, changed] = update_cells_two_signals_multiply_finite_Hill(cells, distances, M_int, a0,...
            Rcell, Con, Coff, K, lambda, hill, noise);

    while t<gz+1
        t = t+1;
        cells = cellsOut;
        cells_hist{end+1} = cells; %{cells(:, 1), cells(:, 2)};
        [period, t_onset] = periodicity_test_short(cells_hist); 
        if display_fig
            pause(0.01);
            update_figure_periodic_scatter(h_cells, cells, t, disp_mol, showI, a0, distances);
            %update_cell_figure_continuum(hin, pos, cells_in, cell_type, t, disp_mol);
        end
        [cellsOut, changed] = update_cells_two_signals_multiply_finite_Hill(cells, distances, M_int, a0,...
            Rcell, Con, Coff, K, lambda, hill, noise);
    end
    
    t_out = t; %default t_out
    
    % if periodicity found, refine check to find period
    if period<Inf
        t_check = 0;
        [period, t_onset] = periodicity_test_detailed(cells_hist, t_check, period);
        t_out = t_onset + period;
    end
    fprintf('Final: t_out = %d, period %d \n', t_out, period);
    %----------------------------------------------------------------------
    % Save result
    if save_folder % Check whether to save the result
        I_ini_str = '';
        if InitiateI
            I_ini_str = sprintf('_I_ini_%.2f_%.2f', I0(1), I0(2));
        end

        fname_str = strrep(sprintf('%s_t_out_%d_period_%d', fname_str_template,...
            t_out, period), '.', 'p');
        ext = '.mat';

        % check if filename already exists
        i = 1;
        fname = fullfile(save_folder, strcat(fname_str, '-v', num2str(i), ext));
        while exist(fname, 'file') == 2
            i=i+1;
            fname = fullfile(save_folder, strcat(fname_str, '-v', num2str(i), ext));
        end
        % variables to save
        rcell = Rcell/a0;      
        lambda12 = lambda(2);
        save_vars = {N, a0, K, Con, Coff, M_int, hill,...
            noise, p0, I0, rcell, Rcell, lambda12, lambda,...
            sim_ID, I_ini_str, tmax, mcsteps};
        save_vars_lbl = {'N', 'a0', 'K', 'Con', 'Coff', 'M_int', 'hill',...
            'noise', 'p0', 'I0', 'rcell', 'Rcell',  'lambda12', 'lambda', ...
            'sim_ID', 'I_ini_str', 'tmax', 'mcsteps'};

        if InitiateI
            save_vars{end+1} = I0;
            save_vars_lbl{end+1} = 'I_ini';
        end

        save_consts_struct = cell2struct(save_vars, save_vars_lbl, 2);

        % save
        save(fname, 'save_consts_struct', 'cells_hist', 't_out',...
            'changed', 'period', 't_onset', 'positions', 'distances');
        fprintf('Saved simulation: %s ; \n', fname);
    end
end