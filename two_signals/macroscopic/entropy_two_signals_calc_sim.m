function [Omega_E, S, frac_E] = entropy_two_signals_calc_sim(N, a0, rcell,...
    lambda, hill, noise, M_int, Con, Coff, K, dist, parent_folder)
    % Check whether simulation has already been done
    % If done, load saved data
    M_int_s = sprintf('%d_%d_%d_%d', M_int(1,1), M_int(1,2), M_int(2,1), M_int(2,2));
    a0_s = sprintf('%.2f', a0);
    R_s = sprintf('%.2f', rcell);
    K_s = sprintf('%d_%d_%d_%d', K(1,1), K(1,2), K(2,1), K(2,2));
    Con_s = sprintf('%d_%d', Con(1), Con(2));
    lambda_s = sprintf('%.1f', lambda(2));
    
    fname_str_base = strrep(sprintf(...
        'N%d_M_int_%s_a0_%s_rcell_%s_K_%s_Con_%s_hill%.1f_noise%.1f',...
         N, M_int_s, a0_s, R_s, K_s, Con_s, hill, noise), '.', 'p');
    
    % store batch data in:
    subfolder1 = strrep(sprintf(...
        'N%d_M_int_%s_a0_%s_rcell_%s_lambda12_%s_Con_%s_hill%.1f_noise%.1f',...
         N, M_int_s, a0_s, R_s, lambda_s, Con_s, hill, noise), '.', 'p'); 
    batch_file_name = strcat('K_', K_s);
    % store individual data files in:
    subfolder2 = strcat('K_', K_s);
     
    if exist(fullfile(parent_folder, subfolder1, ...
            strcat(batch_file_name, '_all.mat')), 'file')==2
        load(fullfile(parent_folder, subfolder1, ...
            strcat(batch_file_name, '_all.mat')), 'Omega_E', 'S', 'frac_E');
        fprintf('Loaded %s \n', fname_str_base');
        %pause(0.1);
        return;
    end 
    %%
    % If not done, proceed with simulation
    % Count total # simulations
    count_sets = nchoosek(N+4-1, N);
    
    % Simulation
    Omega_E_all = [];
    counter = 0;
    for n00=0:N
        %disp(n00);
        for n01=0:N-n00
            for n10=0:N-n00-n01
                n = [n00 n01; n10 N-n00-n01-n10];
                
                counter = counter + 1;
                fprintf('Total # initial conditions: %d, counter = %d \n', count_sets, counter);
                
                % file name
                iniON_s = sprintf('%d_%d_%d_%d', n(1,1), n(1,2), n(2,1), n(2,2));
                fname_str = strrep(sprintf(...
                    '%s_iniON_%s', fname_str_base, iniON_s), '.', 'p');  
                %disp(fname_str);
                
                % skip done files
                folder = fullfile(parent_folder, subfolder1, subfolder2);
                if exist(fullfile(folder, strcat(fname_str, '.mat')), 'file')==2
                    continue;
                end
                
                % multinomial coefficient
                mult_coeff = prod([1:N]./[1:n(1) 1:n(2) 1:n(3) 1:n(4)]); 
                % mult_coeff_2 = nchoosek(n(1)+n(2), n(2))*nchoosek(n(1)+n(2)+n(3),
                % n(3))*nchoosek(n(1)+n(2)+n(3)+n(4), n(4)); % alternative
                
                % Calculate Peq by simulations
                max_sim = 10^4;
                Q = 0; % number of equilibrium states
                nsim = round(min(max_sim, mult_coeff));
                fprintf('Simulations = %d \n', nsim);
                parfor i=1:nsim
                    %disp(i);
                    % initiate random config
                        % subdivide indices for different cell states
                    idx = cell(3,1);
                    idx_not00 = randperm(N, N-n(1,1)); % cell indices of all cells not (0,0)
                    part_idx = [1 n(1,2) n(1,2)+n(2,1) n(1,2)+n(2,1)+n(2,2)]; % partition indices
                    idx{1} = idx_not00(part_idx(1):part_idx(2));
                    idx{2} = idx_not00(part_idx(2)+1:part_idx(3));
                    idx{3} = idx_not00(part_idx(3)+1:part_idx(4));

                    %{
                    %Check that sets have the right number of elems
                    check = (numel(idx{1})==iniON(1,2));
                    check = check*(numel(idx{2})==iniON(2,1));
                    check = check*(numel(idx{3})==iniON(2,2));
                    % Check that sets have disjoint elements
                    check = check*(numel(union(union(idx{1}, idx{2}), idx{3}))==N-iniON(1,1)); 
                    fprintf('Check passed? %d \n', check);
                    %}

                    % turn genes ON
                    cells = zeros(N, 2);
                    cells(idx{1}, 2) = 1;
                    cells(idx{2}, 1) = 1;
                    cells(idx{3}, 1) = 1;
                    cells(idx{3}, 2) = 1;

                    % check if stable
                    Rcell = rcell*a0;
                    [~, changed] = update_cells_two_signals_multiply_finite_Hill(...
                        cells, dist, M_int, a0, Rcell, Con, Coff, K, lambda, hill, noise);
                    if ~changed
                        Q=Q+1;
                    end
                end
                Peq = Q/nsim;

                % Compute Omega_E
                Omega_En = mult_coeff*Peq;
                Omega_E_all(end+1) = Omega_En;
                
                % Save partial results
                save(fullfile(parent_folder, subfolder1, subfolder2, fname_str), 'Omega_En', 'Peq', 'nsim');
            end
        end
    end

    Omega_E = sum(Omega_E_all);
    fprintf('Omega_E = %.3f \n', Omega_E);
    %% Get fraction of equilibrium states
    S = log(Omega_E);
    Smax = N*log(2);
    frac_E = exp(S - Smax);
    fprintf('Fraction equilibrium = %.3f \n', frac_E);
    
    % Save combined results
    fname_str = strcat('K_', K_s, '_all');
    save(fullfile(parent_folder, subfolder1, fname_str),...
        'Omega_E', 'S', 'frac_E', 'N', 'a0', 'rcell',...
        'lambda', 'hill', 'noise', 'M_int', 'Con', 'Coff', 'K');
    %% Check all Omega_E
    % figure;
    % histogram(Omega_E_all);
end
