function [dH_count, dI, dH_mean, dH_std] = ...
    spin_flip_dist_eq_parallel_fixed_pin_hill(dist, Con, K, a0, Rcell, p0, hill, prec, noise)
% generate configurations that differ from random initial configuration by
% 'flip' spins and measure the distance in final configuration
% variable flip: number of cells to flip at each time point

N = size(dist,1); % number of cells
n_smpl = 100;
nflips = 100; % number of deviating lattices per random initial lattice
dH_count = zeros(N+1,1);

if p0 == 0 || p0==1 %special case all cells ON/OFF, only 1 base configuration
    dI = zeros(1, nflips);
    dH = zeros(1, nflips);
    % Base case
    cells = p0*ones(N,1); 
    [cells, changed] = update_cells_continuum(cells, dist, Con, K, a0, Rcell, hill, prec);
    while changed
        [cells, changed] = update_cells_continuum(cells, dist, Con, K, a0, Rcell, hill, prec);
    end
    cells_base_out = cells;

    % Run lattice with initial noise
    for idx = 1:nflips 
        cells = p0*ones(N,1);
        cells = generate_pattern_noise(cells, noise); % add noise to original pattern 
        [cells, changed] = update_cells_continuum(cells, dist, Con, K, a0, Rcell, hill, prec);
        while changed
            [cells, changed] = update_cells_continuum(cells, dist, Con, K, a0, Rcell, hill, prec);
        end
        % compute Hamming distance
        dH_out = sum(abs(cells-cells_base_out));
        dH(idx) = dH_out;
        dH_out_idx = round(dH_out);
        dH_count(dH_out_idx+1) = dH_count(dH_out_idx+1)+1;
    end
    dH_mean = mean(dH);
    dH_std = std(dH);
else %other cases with intermediate number of ON cells
    dI = zeros(n_smpl, nflips);
    dH_aux = zeros(n_smpl, nflips);
    % Test n_smpl random initial configurations and update cells
    % until they are in final configuration
    parfor i = 1:n_smpl
        disp(i);
        % base case
        cells = init_random_cells_montecarlo(N, p0);
        cells_base_ini = cells; % store base initial config
        [cells, changed] = update_cells_continuum(cells, dist, Con, K, a0, Rcell, hill, prec);
        while changed
            [cells, changed] = update_cells_continuum(cells, dist, Con, K, a0, Rcell, hill, prec);
        end
        cells_base_out = cells;
        [I_base, ~] = moranI(cells, a0*dist);

        % Run lattice with initial noise
        for idx = 1:nflips 
            cells = generate_pattern_noise(cells_base_ini, noise); % add noise to original pattern
            [cells, changed] = update_cells_continuum(cells, dist, Con, K, a0, Rcell, hill, prec);
            while changed
                [cells, changed] = update_cells_continuum(cells, dist, Con, K, a0, Rcell, hill, prec);
            end
            % compute distance and dI
            dH_aux(i, idx) = sum(abs(cells-cells_base_out));
            [I_out, ~] = moranI(cells, a0*dist); % final I
            dI(i, idx) = I_out - I_base;
        end
    end    
    % store result into dH_count
    dH_aux_round = round(dH_aux);
    for dH_out = 0:N
        dH_count(dH_out+1) = dH_count(dH_out+1) + sum(sum(dH_aux_round == dH_out));
    end
    dH_mean = mean(mean(dH_aux));
    std_temp = std(dH_aux, 0, 2);
    dH_std = sqrt(sum(std_temp.^2))/n_smpl;
end

end
