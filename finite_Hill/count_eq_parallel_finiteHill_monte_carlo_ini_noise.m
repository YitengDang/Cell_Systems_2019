function [peqcount, dHcount, t_av, I_av, dH_mean, dH_std] =...
    count_eq_parallel_finiteHill_monte_carlo_ini_noise(dist, Con, K, a0, Rcell, hill, prec, noise)
% Count the number of equilibrium states found in n_smpl random
% configurations trials, using parallel processing, making the simulation
% faster for large systems.

% This is used to simulate the p_ini vs. p_eq map. It starts with random
% configurations and runs until equilibrium, saving a 2d map, count, that
% for each p_ini and p_eq counts the number of simulations that fall within
% it. t_av counts the average number of steps that it takes to reach
% equilibrium.
N = size(dist,1); % number of cells
n_smpl = 1000; % # random initial configs
noise_smpl = 1; % #samples with added noise
peqcount = zeros(N+1,N+1);
dHcount = zeros(N+1,N+1);
t_av = zeros(N+1,1);
I_av = zeros(N+1,3);
dH_mean = zeros(N+1, 1);
dH_std = zeros(N+1, 1);

for k = 0 : N
    disp(k)
    if k == 0 || k==N %special case all cells OFF
        % Base case
        cells = k/N*ones(N,1); 
        [cells, changed] = update_cells_continuum(cells, dist, Con, K, a0, Rcell, hill, prec);
        while changed
            t_av(k+1) = t_av(k+1) + 1;
            [cells, changed] = update_cells_continuum(cells, dist, Con, K, a0, Rcell, hill, prec);
        end
        cells_base_out = cells;
        kout = round(sum(cells));
        peqcount(k+1,kout+1) = peqcount(k+1,kout+1) + 1;
        
        % Spin flip (N possibilities)
        dH_temp = zeros(N, 1);
        for flip_idx = 1:n_smpl*noise_smpl
            cells = k/N*ones(N,1);
            cells = generate_pattern_noise(cells, noise); % add noise to original pattern
            [cells, changed] = update_cells_continuum(cells, dist, Con, K, a0, Rcell, hill, prec);
            while changed
                [cells, changed] = update_cells_continuum(cells, dist, Con, K, a0, Rcell, hill, prec);
            end
            % compute Hamming distance
            dH_out = round(sum(abs(cells-cells_base_out)));
            dH_temp(flip_idx) = dH_out;
            dHcount(k+1,dH_out+1) = dHcount(k+1,dH_out+1)+1;
        end
        
        dH_mean(k+1) = mean(dH_temp);
        dH_std(k+1) = std(dH_temp);
    else %other cases with intermediate number of ON cells
        t_aux = zeros(n_smpl,1);
        I_aux = zeros(n_smpl,1);
        k_aux = zeros(n_smpl,1);
        dH_aux = zeros(n_smpl,N);
        % Test n_smpl random initial configurations and update cells
        % until they are in final configuration
        parfor i = 1:n_smpl
            cells_base_ini = init_random_cells_montecarlo(N, k/N);
            [cells, changed] = update_cells_continuum(cells_base_ini, dist, Con, K, a0, Rcell, hill, prec);
            while changed
                t_aux(i) = t_aux(i) + 1;
                [cells, changed] = update_cells_continuum(cells, dist, Con, K, a0, Rcell, hill, prec);
            end
            cells_base_out = cells;
            kout = round(sum(cells));
            k_aux(i) = kout;
            [I_aux(i), ~] = moranI(cells, a0*dist);
            
            % Spin flip (N possibilities)
            for flip_idx = 1:noise_smpl
                cells = cells_base_ini;
                cells = generate_pattern_noise(cells, noise); % add noise to original pattern
                [cells, changed] = update_cells_continuum(cells, dist, Con, K, a0, Rcell, hill, prec);
                while changed
                    [cells, changed] = update_cells_continuum(cells, dist, Con, K, a0, Rcell, hill, prec);
                end
                % compute Hamming distance
                dH_aux(i, flip_idx) = sum(abs(cells-cells_base_out));
            end
        end
        dH_mean(k+1) = mean(mean(dH_aux));
        std_temp = std(dH_aux, 0, 2);
        dH_std(k+1) = sqrt(sum(std_temp.^2))/N;
        t_av(k+1) = sum(sum(t_aux))/n_smpl;
        for i = 0:N
            peqcount(k+1,i+1) = peqcount(k+1,i+1) + sum(sum(k_aux == i));
            dHcount(k+1,i+1) = dHcount(k+1,i+1) + sum(sum(round(dH_aux) == i));
        end
        I_av(k+1,1) = sum(sum(I_aux))/n_smpl;
        I_av(k+1,2) = min(I_aux);
        I_av(k+1,3) = max(I_aux);
    end
end

end
