function [cells_hist, period, t_onset] = time_evolution_save_func_efficient_checks_change_cells(...
    N, a0, Rcell, lambda, hill, noise, M_int, K, Con, Coff,...
    distances, positions, sim_ID, mcsteps, tmax,...
    fname_save, cells_ini, num_cells_changed)
    % Similar to time_evolution_save_func, but
    % implements an efficient periodicity test by checking only every
    % t_check time steps whether trajectory has become periodic (after a
    % transient phase of t_ac steps when it checks every step).
    
    % Starts from an input initial configuration and changes a number of
    % cells
    % cells_ini: initial config
    % num_flip_cells: number of cells to flip
    %%
    % Simulate
    cells_hist = {};
    
    % store initial config
    rand_idx = randperm(N, num_cells_changed); % choose random cell
    switch randi(3)
        case 1 % flip state 1
            cells_ini(rand_idx, 1) = ~cells_ini(rand_idx, 1);
        case 2 % flip state 2
            cells_ini(rand_idx, 2) = ~cells_ini(rand_idx, 2);
        case 3 % flip both states
            cells_ini(rand_idx, :) = ~cells_ini(rand_idx, :);
    end
    cells_hist{end+1} = cells_ini; %{cells(:, 1), cells(:, 2)};
    
    p0 = mean(cells_ini, 1);
    I0(1) = moranI(cells_ini(:,1), a0*distances);
    I0(2) = moranI(cells_ini(:,2), a0*distances);
    
    %-------------dynamics-----------------------------------------
    t = 0;
    period = Inf; %default values
    t_onset = Inf; 
    [cellsOut, changed] = update_cells_two_signals_multiply_finite_Hill(cells_ini, distances, M_int, a0,...
            Rcell, Con, Coff, K, lambda, hill, noise);

    t_ac = 10^2; 
    while changed && period==Inf && t<t_ac
        t = t+1;
        cells = cellsOut;
        cells_hist{end+1} = cells; %{cells(:, 1), cells(:, 2)};
        [period, t_onset] = periodicity_test_short(cells_hist); 
        %update_cell_figure_continuum(app, pos, dist, a0, cells, app.Time, cell_type, disp_mol, 0);
        [cellsOut, changed] = update_cells_two_signals_multiply_finite_Hill(cells, distances, M_int, a0,...
            Rcell, Con, Coff, K, lambda, hill, noise);
    end
    
    % check periodically after t_ac time steps, with period t_check
    t_check = 10^3; 
    while changed && period==Inf && t<tmax
        t = t+1;
        cells = cellsOut;
        cells_hist{end+1} = cells; %{cells(:, 1), cells(:, 2)};
        if mod(t, t_check)==0
            [period, t_onset] = periodicity_test_short(cells_hist); 
        end
        %update_cell_figure_continuum(app, pos, dist, a0, cells, app.Time, cell_type, disp_mol, 0);
        [cellsOut, changed] = update_cells_two_signals_multiply_finite_Hill(cells, distances, M_int, a0,...
            Rcell, Con, Coff, K, lambda, hill, noise);
    end

    t_out = t; %default t_out
    % if periodicity found, refine check to find period
    if period<Inf && t>t_ac
        [period, t_onset] = periodicity_test_detailed(cells_hist, t_check, period);
        t_out = t_onset + period;
    end

    %{
    tmax_string = '';
    if changed && t==tmax
        tmax_string = '_tmax_reached';
    end
    %}

    fprintf('Final: t_out = %d, period %d \n', t_out, period);
    %----------------------------------------------------------------------
    % Save result
    %fname_str = strrep(sprintf('two_signal_mult_N%d_K12_%d_t_out_%d_period_%d_%dcells_changed',...
    %	N, K(1,2), t_out, period, num_cells_changed), '.', 'p');
    
    %{        
    labels = {'_tmax_reached', ''};
    if t_out==tmax
        label_choice = 1;
    else
        label_choice = 2;
    end
    %}
    ext = '.mat';
    
    % check if filename already exists
    i=1;
    fname = strcat(fname_save, '-v', num2str(i), ext);
    while exist(fname, 'file') == 2
        i=i+1;
        fname =  strcat(fname_save, '-v', num2str(i), ext);
    end
    
    % variables to save
    I_ini_str = '';
    InitiateI = 0;
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