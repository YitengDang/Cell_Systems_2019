function [dH_count, dI, dH_mean, dH_std] = ...
    spin_flip_dist_eq_parallel_fixed_pin(dist, Con, K, a0, Rcell, k, flip)
% generate configurations that differ from random initial configuration by
% 'flip' spins and measure the distance in final configuration
% variable flip: number of cells to flip at each time point

N = size(dist,1); % number of cells
n_smpl = 100;
nflips = 1000; % number of random flips per lattice
dH_count = zeros(N+1,1);

if k == 0 || k==N %special case all cells ON/OFF, only 1 base configuration
    dI = zeros(1, nflips);
    dH = zeros(1, nflips);
    % Base case
    cells = k/N*ones(N,1); 
    [cells, changed] = update_cells(cells, dist, Con, K, a0, Rcell);
    while changed
        [cells, changed] = update_cells(cells, dist, Con, K, a0, Rcell);
    end
    cells_base_out = cells;

    % Spin flip (N possibilities)
    for idx = 1:nflips 
        cells = k/N*ones(N,1);
        idx_flip = randperm(N, flip);
        cells(idx_flip) = ~cells(idx_flip);
        [cells, changed] = update_cells(cells, dist, Con, K, a0, Rcell);
        while changed
            [cells, changed] = update_cells(cells, dist, Con, K, a0, Rcell);
        end
        % compute Hamming distance
        dH_out = sum(abs(cells-cells_base_out));
        dH(idx) = dH_out;
        dH_count(dH_out+1) = dH_count(dH_out+1)+1;
    end
    dH_mean = mean(dH);
    dH_std = std(dH);
else %other cases with intermediate number of ON cells
    dI = zeros(n_smpl, nflips);
    dH_aux = zeros(n_smpl, nflips);
    % Test n_smpl random initial configurations and update cells
    % until they are in final configuration
    parfor i = 1:n_smpl
        % base case
        rnd_idx = randperm(N,k);
        cells = zeros(N,1);
        cells(rnd_idx) = 1;
        cells_base_ini = cells; % store base initial config
        [cells, changed] = update_cells(cells, dist, Con, K, a0, Rcell);
        while changed
            [cells, changed] = update_cells(cells, dist, Con, K, a0, Rcell);
        end
        cells_base_out = cells;
        [I_base, ~] = moranI(cells, a0*dist);

        % Spin flip (N possibilities)
        for idx = 1:nflips 
            cells = cells_base_ini; 
            idx_flip = randperm(N, flip);
            cells(idx_flip) = ~cells(idx_flip);
            [cells, changed] = update_cells(cells, dist, Con, K, a0, Rcell);
            while changed
                [cells, changed] = update_cells(cells, dist, Con, K, a0, Rcell);
            end
            % compute distance and dI
            dH_aux(i, idx) = sum(abs(cells-cells_base_out));
            [I_out, ~] = moranI(cells, a0*dist); % final I
            dI(i, idx) = I_out - I_base;
        end
    end    
    % store result into dH_count
    for dH_out = 0:N
        dH_count(dH_out+1) = dH_count(dH_out+1) + sum(sum(dH_aux == dH_out));
    end
    dH_mean = mean(mean(dH_aux));
    std_temp = std(dH_aux, 0, 2);
    dH_std = sqrt(sum(std_temp.^2))/n_smpl;
end

end
