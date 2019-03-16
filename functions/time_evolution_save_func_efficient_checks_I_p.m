function [cells_hist, period, t_onset] = time_evolution_save_func_efficient_checks_I_p(...
    N, a0, Rcell, lambda, hill, noise, M_int, K, Con, Coff,...
    distances, positions, sim_ID, mcsteps, InitiateI, p0, I0, tmax, save_folder,sim_it)
    % Similar to time_evolution_save_func, but
    % implements an efficient periodicity test by checking only every
    % t_check time steps whether trajectory has become periodic (after a
    % transient phase of t_ac steps when it checks every step).

    % Simulate
    cells_hist = {};
    
    % generate initial lattice
    if p0==Inf % special option: generate random lattice
        cells = randi(2, N, 2)-1; 
    else
        iniON = round(p0*N);
        cells = zeros(N, 2);
        for i=1:numel(iniON)
            cells(randperm(N,iniON(i)), i) = 1;
            if InitiateI && hill==Inf
                %fprintf('Generating lattice with I%d(t=0)... \n', i);
                dI = 0.05;
                [cells_temp, ~, ~] = generate_I_new(cells(:, i), I0(i), I0(i)+dI, distances, a0,Rcell);
                cells(:,i) = cells_temp;
                %fprintf('Generated initial I%d: %.2f; in range (Y=1/N=0)? %d; \n', i, I_ini, test);
            end
        end
    end
    % store initial config
    cells_hist{end+1} = cells; %{cells(:, 1), cells(:, 2)};
    %-------------dynamics-----------------------------------------
    t = 0;
    period = Inf; %default values
    t_onset = Inf; 
    [cellsOut, changed] = update_cells_two_signals_multiply_finite_Hill(cells, distances, M_int, a0,...
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

    tmax_string = '';
    if changed && t==tmax
        tmax_string = '_tmax_reached';
    end

    fprintf('Final: t_out = %d, period %d \n', t_out, period);
    %----------------------------------------------------------------------
    % Save result
    
    I_ini_str = '';
    if InitiateI
        I_ini_str = sprintf('I_ini_%.2f_%.2f', I0(1), I0(2));
    end
    
    p_ini_str = sprintf('p_ini_%.2f_%.2f', p0(1), p0(2));
    
    fname_str = strrep(sprintf('%s_sim_%g_%s_%s_t_out_%g_period_%s',...
                sim_ID,sim_it, I_ini_str,p_ini_str, t_out, num2str(period)), '.', 'p');
    %fname_str = strrep(sprintf(...
    %    'N%d_initiateI%d%s_K12_%.2f_K21_%.2f_Con1_%.2f_Con2_%.2f_t_out_%d_period_%s',...
    %    N, InitiateI, I_ini_str, K(1,2), K(2,1), Con(1), Con(2),...
    %    t_out, num2str(period)), '.', 'p');
    %fname_str = 'all_topologies_simulate';
    
    labels = {'_tmax_reached', ''};
    if t_out==tmax
        label_choice = 1;
    else
        label_choice = 2;
    end
    ext = '.mat';
    
    % check if filename already exists
    i=1;
    fname_all = cell(2,1); % all possible names
    fname_all{1} = fullfile(save_folder, strcat(fname_str, '-v', num2str(i), labels{1}, ext));
    fname_all{2} = fullfile(save_folder, strcat(fname_str, '-v', num2str(i), labels{2}, ext)); % check both labels 
    while exist(fname_all{1}, 'file') == 2 || exist(fname_all{2}, 'file') == 2
        i=i+1;
        fname_all{1} = fullfile(save_folder, strcat(fname_str, '-v', num2str(i), labels{1}, ext));
        fname_all{2} = fullfile(save_folder, strcat(fname_str, '-v', num2str(i), labels{2}, ext));
    end
    fname = fname_all{label_choice};
    
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