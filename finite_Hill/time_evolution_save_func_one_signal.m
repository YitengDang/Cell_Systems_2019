function [cells_hist, period, t_onset] = time_evolution_save_func_one_signal(...
    N, a0, Rcell, hill, noise, K, Con, prec,...
    distances, positions, sim_ID, mcsteps, InitiateI, p_ini, I_ini, cells_ini, tmax,...
    save_folder, fname_str_template, display_fig)
    % Similar to time_evolution_save_func, but
    % implements an efficient periodicity test by checking only every
    % t_check time steps whether trajectory has become periodic (after a
    % transient phase of t_ac steps when it checks every step).
    cells_hist = {};
    gz = sqrt(N);
    rcell = Rcell/a0;
    
    % randomize initial lattice
    if isempty(positions) && isempty(distances)
        nodisplay = 1;
        [positions, distances] = initial_cells_random_markov_periodic(gz, mcsteps,...
            rcell, nodisplay);
    end
    
    % generate initial lattice
    if ~isempty(cells_ini)
        cells = cells_ini;
    else
        if p_ini==Inf % special option: generate random lattice
            if hill==Inf
                cells = randi(2, N, 1)-1; % binary cells
            else
                % cells = rand(N, 2); % continuous cells, uniformly sampled
                sigma = 1/4;
                cells = (randn(N, 1))*sigma + 1/2; % continuous cells, normal distribution
                cells(cells>1) = 1;
                cells(cells<0) = 0;           
            end
        else
            iniON = round(p_ini*N);
            cells = zeros(N, 1);
            cells(randperm(N, iniON)) = 1;
            if InitiateI && hill==Inf
                %fprintf('Generating lattice with I%d(t=0)... \n', i);
                dI = 0.1;
                [cells, ~, ~] = generate_I_new(cells, I_ini, I_ini+dI, distances, a0);
                %fprintf('Generated initial I%d: %.2f; in range (Y=1/N=0)? %d; \n', i, I_ini, test);
            end
        end
    end
    % store initial config
    cells_hist{end+1} = cells; %{cells(:, 1), cells(:, 2)};
    
    % Show initial condiguration
    if display_fig
        hin = figure;
        %set(hin, 'Position', [100 100 800 800]);
        [h_cells, h_borders] = reset_cell_figure(hin, positions, rcell);
        disp_mol = 1;
        showI = 0; 
        t=0;
        %disp(size(cells));
        update_figure_periodic_scatter(h_cells, cells, t, disp_mol, showI, a0, distances);
        pause(1);
    end
    
    %-------------dynamics-------------------------------------------------
    t = 0;
    [cells_out, changed, ~] = ...
            update_cells_continuum(cells, distances, Con, K, a0, Rcell, hill, prec);
    while changed && t<tmax
        t = t+1;
        cells = cells_out;
        cells_hist{end+1} = cells; %{cells(:, 1), cells(:, 2)};
        if display_fig
            pause(0.01);
            update_figure_periodic_scatter(h_cells, cells, t, disp_mol, showI, a0, distances);
            %update_cell_figure_continuum(hin, pos, cells_in, cell_type, t, disp_mol);
        end
        [cells_out, changed, ~] = ...
            update_cells_continuum(cells, distances, Con, K, a0, Rcell, hill, prec);
    end
    t_out = t; %default t_out
    fprintf('Final: t_out = %d \n', t_out);
    %----------------------------------------------------------------------
    % Save result
    I_ini_str = '';
    if InitiateI
        I_ini_str = sprintf('_I_ini_%.2f', I_ini);
    end
    
    % Format of output filename
    fname_str = strrep(sprintf('%s_t_out_%d', fname_str_template,...
        t_out), '.', 'p');
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
    M_int = 1;
    Coff = 1;
    p_ini = mean(cells_hist{1});
    I_ini = moranI(cells_hist{1}, a0*distances);
    lambda = 1;
    save_vars = {N, a0, K, Con, Coff, M_int, hill,...
        noise, p_ini, I_ini, rcell, Rcell, lambda,...
        sim_ID, I_ini_str, tmax, mcsteps};
    save_vars_lbl = {'N', 'a0', 'K', 'Con', 'Coff', 'M_int', 'hill',...
        'noise', 'p_ini', 'I_ini', 'rcell', 'Rcell', 'lambda', ...
        'sim_ID', 'I_ini_str', 'tmax', 'mcsteps'};
    
    if InitiateI
        save_vars{end+1} = I_ini;
        save_vars_lbl{end+1} = 'I_ini';
    end

    save_consts_struct = cell2struct(save_vars, save_vars_lbl, 2);
    
    % save
    save(fname, 'save_consts_struct', 'cells_hist', 't_out',...
        'changed', 'positions', 'distances');
    fprintf('Saved simulation: %s ; \n', fname);
end