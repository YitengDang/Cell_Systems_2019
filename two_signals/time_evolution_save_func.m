function [cells_hist, period, t_onset] = time_evolution_save_func(...
    N, a0, Rcell, lambda, hill, noise, M_int, K, Con, Coff,...
    dist, p0, InitiateI, save_folder)

    %% Simulate
    cells_hist = {};

    % generate initial lattice
    iniON = round(p0*N);
    cells = zeros(N, 2);
    for i=1:numel(iniON)
        cells(randperm(N,iniON(i)), i) = 1;
        if InitiateI && hill==Inf
            %fprintf('Generating lattice with I%d(t=0)... \n', i);
            dI = 0.1;
            [cells_temp, ~, ~] = generate_I_new(cells(:, i), I0(i), I0(i)+dI, dist, a0);
            cells(:,i) = cells_temp;
            %fprintf('Generated initial I%d: %.2f; in range (Y=1/N=0)? %d; \n', i, I_ini, test);
        end
    end

    % store initial config
    cells_hist{end+1} = cells; %{cells(:, 1), cells(:, 2)};
    %-------------dynamics-----------------------------------------
    t = 0;
    period = Inf; %default values
    t_onset = Inf; 
    [cellsOut, changed] = update_cells_two_signals_multiply_finite_Hill(cells, dist, M_int, a0,...
            Rcell, Con, Coff, K, lambda, hill, noise);

    while changed && period==Inf %&& t<tmax
        t = t+1;
        cells = cellsOut;
        cells_hist{end+1} = cells; %{cells(:, 1), cells(:, 2)};
        [period, t_onset] = periodicity_test_short(cells_hist); 
        [cellsOut, changed] = update_cells_two_signals_multiply_finite_Hill(cells, dist, M_int, a0,...
            Rcell, Con, Coff, K, lambda, hill, noise);
    end
    t_out = t; % save final time
    fprintf('t_out = %d, period %d \n', t_out, period);
    
    %----------------------------------------------------------------------
    % Save result
    
    I_ini_str = '';
    if InitiateI
        I_ini_str = sprintf('_I_ini_%.2f_%.2f', I0(1), I0(2));
    end
        
    %fname_str = strrep(sprintf(...
    %    'N%d_initiateI%d%s_K12_%.2f_K21_%.2f_Con1_%.2f_Con2_%.2f_t_out_%d_period_%s',...
    %    N, InitiateI, I_ini_str, K(1,2), K(2,1), Con(1), Con(2),...
    %    t_out, num2str(period)), '.', 'p');
    fname_str = 'all_topologies_simulate';
    
    ext = '.mat';
    label = '';

    % check if filename already exists
    i=1;
    fname = fullfile(save_folder, strcat(fname_str, '-v', num2str(i), label, ext));
    while exist(fname, 'file') == 2
        i=i+1;
        fname = fullfile(save_folder, strcat(fname_str, '-v', num2str(i), label, ext));
    end

    if InitiateI
        save_vars = {N, a0, K, Con, Coff, M_int, hill, noise, p0, I0, Rcell,...
            lambda, I_ini_str};
        save_vars_lbl = {'N', 'a0', 'K', 'Con', 'Coff', 'M_int', 'hill',...
            'noise', 'p_ini', 'I_ini', 'Rcell',...
            'lambda', 'I_ini_str'};
    else
        save_vars = {N, a0, K, Con, Coff, M_int, hill, noise, p0, Rcell,...
            lambda, I_ini_str};
        save_vars_lbl = {'N', 'a0', 'K', 'Con', 'Coff', 'M_int', 'hill',...
            'noise', 'p_ini', 'Rcell',...
            'lambda', 'I_ini_str'};
    end

    save_consts_struct = cell2struct(save_vars, save_vars_lbl, 2);

    save(fname, 'save_consts_struct', 'cells_hist', 't_out',...
        'changed', 'period', 't_onset');
        fprintf('Saved simulation: %s ; \n', fname);
end