function [count, t_av, I_av, ntmax] = count_eq_parallel_finiteHill_monte_carlo(...
    dist, Con, K, Rcell, a0, hill, prec)
% Count the number of equilibrium states found in n_smpl random
% configurations trials, using parallel processing, making the simulation
% faster for large systems.

% This is used to simulate the p_ini vs. p_eq map. It starts with random
% configurations and runs until equilibrium, saving a 2d map, count, that
% for each p_ini and p_eq counts the number of simulations that fall within
% it. t_av counts the average number of steps that it takes to reach
% equilibrium.

% k runs from 0 to N, <Xi>(initial) = k/N
% Uniform sampling: start with a lattice of identical cells (Xi = p for all
% i). Only one simulation needed without noise.

N = size(dist,1); % number of cells
n_smpl = 1000; % number of samples for each value
count = zeros(N+1,N+1); % counts of #eq. states
t_av = zeros(N+1,1);
I_av = zeros(N+1,3);
tmax = 2000;
ntmax = zeros(N+1, 1); % number of simulations reaching tmax
for k = 0:N % k = <Xi(t=0)> = p_initial
    disp(k)
    if k == 0 %special case all cells OFF; only one state -> one simulation enough
        cells = zeros(N,1);
        [cells, changed] = update_cells_continuum(cells, dist, Con, K, a0, Rcell, hill, prec);
        while changed && t_av(k+1)<=tmax
            t_av(k+1) = t_av(k+1) + 1;
            [cells, changed] = update_cells_continuum(cells, dist, Con, K, a0, Rcell, hill, prec);
        end
        %--- count tmax---
        if t_av(k+1)==tmax
            ntmax(k+1) = ntmax(k+1);
        end
        %-----------------
        kout = round(sum(cells));
        I_av(k+1, :) = repmat(moranI(cells, a0*dist), 1, 3);
        count(k+1,kout+1) = count(k+1,kout+1) + 1;
    elseif k == N %special case all cells ON
        cells = ones(N,1); % initialize all cells ON
        [cells, changed] = update_cells_continuum(cells, dist, Con, K, a0, Rcell, hill, prec);
        while changed && t_av(k+1)<=tmax
            t_av(k+1) = t_av(k+1) + 1;
            [cells, changed] = update_cells_continuum(cells, dist, Con, K, a0, Rcell, hill, prec);
        end
        %--- count tmax---
        if t_av(k+1)==tmax
            ntmax(k+1) = ntmax(k+1);
        end
        %-----------------
        kout = round(sum(cells));
        I_av(k+1, :) = repmat(moranI(cells, a0*dist), 1, 3);
        count(k+1,kout+1) = count(k+1,kout+1) + 1;
    else %other cases with intermediate number of ON cells
        t_aux = zeros(n_smpl,1);
        I_aux = zeros(n_smpl,1);
        k_aux = zeros(n_smpl,1);
        %--- count tmax---
        ntmax_aux = zeros(n_smpl,1);
        %-----------------
        % Test n_smpl random initial configurations and update cells
        % until they are in final configuration
        parfor i = 1:n_smpl
            %p = k/N; % <Xi> 
            cells = init_random_cells_montecarlo(N, k/N);
            [cells, changed] = update_cells_continuum(cells, dist, Con, K, a0, Rcell, hill, prec);
            while changed && t_aux(i)<=tmax
                t_aux(i) = t_aux(i) + 1;
                [cells, changed] = update_cells_continuum(cells, dist, Con, K, a0, Rcell, hill, prec);
            end
            %--- count tmax---
            if t_aux(i)==tmax
                ntmax_aux(i) = 1;
            end
            %-----------------
            kout = round(sum(cells));
            k_aux(i) = kout;
            [I_aux(i), ~] = moranI(cells, a0*dist);
        end
        %--- count tmax---
        ntmax(k+1) = sum(ntmax_aux);
        %-----------------
        % update count
        for i = 0:N
            count(k+1,i+1) = count(k+1,i+1) + sum(sum(k_aux == i));
        end
        % mean final t
        t_av(k+1) = sum(sum(t_aux))/n_smpl;
        % mean, min and max final I
        I_av(k+1,1) = sum(sum(I_aux))/n_smpl;
        I_av(k+1,2) = min(I_aux);
        I_av(k+1,3) = max(I_aux);
    end
end
    
end

