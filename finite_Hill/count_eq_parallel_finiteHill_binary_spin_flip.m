function [peqcount, dHcount, t_av, I_av, dH_mean, dH_std] = ...
    count_eq_parallel_finiteHill_binary_spin_flip(dist, Con, K, Rcell, a0, hill, prec)
% Count the number of equilibrium states found in n_smpl random
% configurations trials, using parallel processing, making the simulation
% faster for large systems.

% This is used to simulate the p_ini vs. p_eq map. It starts with random
% configurations and runs until equilibrium, saving a 2d map, count, that
% for each p_ini and p_eq counts the number of simulations that fall within
% it. t_av counts the average number of steps that it takes to reach
% equilibrium.

N = size(dist,1); % number of cells
n_smpl = 1000;
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
        kout = sum(cells);
        peqcount(k+1,kout+1) = peqcount(k+1,kout+1) + 1;
        
        % Spin flip (N possibilities)
        dH_temp = zeros(N, 1);
        for cell_idx = 1:N
            cells = k/N*ones(N,1);
            cells(cell_idx) = ~cells(cell_idx);
            [cells, changed] = update_cells_continuum(cells, dist, Con, K, a0, Rcell, hill, prec);
            while changed
                [cells, changed] = update_cells_continuum(cells, dist, Con, K, a0, Rcell, hill, prec);
            end
            % compute Hamming distance
            dH_out = sum(abs(cells-cells_base_out));
            dH_temp(cell_idx) = dH_out;
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
            % base case
            rnd_idx = randperm(N,k);
            cells = zeros(N,1);
            cells(rnd_idx) = 1;
            cells_base_ini = cells; % store base initial config
            [cells, changed] = update_cells_continuum(cells, dist, Con, K, a0, Rcell, hill, prec);
            while changed
                t_aux(i) = t_aux(i) + 1;
                [cells, changed] = update_cells_continuum(cells, dist, Con, K, a0, Rcell, hill, prec);
            end
            cells_base_out = cells;
            kout = sum(cells);
            k_aux(i) = kout;
            [I_aux(i), ~] = moranI(cells, a0*dist);
            
            % Spin flip (N possibilities)
            for cell_idx = 1:N
                cells = cells_base_ini; 
                cells(cell_idx) = ~cells(cell_idx);
                [cells, changed] = update_cells_continuum(cells, dist, Con, K, a0, Rcell, hill, prec);
                while changed
                    [cells, changed] = update_cells_continuum(cells, dist, Con, K, a0, Rcell, hill, prec);
                end
                % compute Hamming distance
                dH_aux(i, cell_idx) = sum(abs(cells-cells_base_out));
            end
            %dH_mean_aux(i) = mean(dH_temp);
            %dH_std_aux(i) = std(dH_temp);
        end
        %dH_mean(k+1) = mean(dH_mean_aux);
        %dH_std(k+1) = sqrt( sum(dH_std_aux.^2) )/N;
        dH_mean(k+1) = mean(mean(dH_aux));
        std_temp = std(dH_aux, 0, 2);
        dH_std(k+1) = sqrt(sum(std_temp.^2))/N;
        
        t_av(k+1) = sum(sum(t_aux))/n_smpl;
        for i = 0:N
            peqcount(k+1,i+1) = peqcount(k+1,i+1) + sum(sum(k_aux == i));
            dHcount(k+1,i+1) = dHcount(k+1,i+1) + sum(sum(dH_aux == i));
        end
        I_av(k+1,1) = sum(sum(I_aux))/n_smpl;
        I_av(k+1,2) = min(I_aux);
        I_av(k+1,3) = max(I_aux);
    end
end

end
