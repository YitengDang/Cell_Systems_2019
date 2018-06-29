function [count, t_av] = count_eq_parallel_markovMC_temp(n, Con, K, a0, rcell, noise, mcsteps, n_smpl)
% version: Markov MC to generate initial lattice

% Count the number of equilibrium states found in n_smpl random
% configurations trials, using parallel processing, making the simulation
% faster for large systems.

% This is used to simulate the p_ini vs. p_eq map. It starts with random
% configurations and runs until equilibrium, saving a 2d map, count, that
% for each p_ini and p_eq counts the number of simulations that fall within
% it. t_av counts the average number of steps that it takes to reach
% equilibrium.

% dimensions
N = n^2;
Lx = 1; % default
Ly = sqrt(3)/2*Lx;
R = rcell*Lx/n; % cell radius

% output vars
%n_smpl = 1000;
count = zeros(N+1,N+1);
t_av = zeros(N+1,1);
%I_av = zeros(N+1,3);
nodisplay = 1;

% loop over Non(ini)
for k = 0 : N
    disp(k)
    if k == 0 %special case all cells OFF
        cells = zeros(N,1);
        [~, dist] = initial_cells_random_markov_periodic(n, mcsteps, rcell, nodisplay);
        

        %---dynamics---
        t=0;
        [cells, changed] = ...
                update_cells_noise(cells, dist, Con, K, a0, rcell*a0, noise);
        while changed
            [cells, changed] = ...
                update_cells_noise(cells, dist, Con, K, a0, rcell*a0, noise);
            t = t+1;
        end
        %---store vars---
        kout = sum(cells);
        count(k+1,kout+1) = count(k+1,kout+1) + 1;
        t_av(k+1) = t;
    elseif k == N %special case all cells ON
        cells = ones(N,1); % initialize all cells ON
        [~, dist] = initial_cells_random_markov_periodic(n, mcsteps, rcell, nodisplay);
        %---dynamics---
        t=0;
        [cells, changed] = ...
                update_cells_noise(cells, dist, Con, K, a0, rcell*a0, noise);
        while changed
            [cells, changed] = ...
                update_cells_noise(cells, dist, Con, K, a0, rcell*a0, noise);
            t = t+1;
        end
        %---store vars---
        kout = sum(cells);
        count(k+1,kout+1) = count(k+1,kout+1) + 1;
        t_av(k+1) = t;
    else %other cases with intermediate number of ON cells
        t_aux = zeros(n_smpl,1);
        %I_aux = zeros(n_smpl,1);
        k_aux = zeros(n_smpl,1);
        
        % Test n_smpl random initial configurations and update cells
        % until they are in final configuration
        parfor i = 1:n_smpl
            % initiate lattice
            [~, dist] = initial_cells_random_markov_periodic(n, mcsteps, rcell, nodisplay);
            rnd_idx = randperm(N,k);
            cells = zeros(N,1);
            cells(rnd_idx) = 1;
            %---dynamics---
            t=0;
            [cells, changed] = ...
            	update_cells_noise(cells, dist, Con, K, a0, rcell*a0, noise);
            while changed
                [cells, changed] = ...
                    update_cells_noise(cells, dist, Con, K, a0, rcell*a0, noise);
                t = t+1;
            end
            %---store vars---
            kout = sum(cells);
            k_aux(i) = kout;
            t_aux(i) = t;
            %I_aux(i) = moranI(cells, a0*dist);
        end
        t_av(k+1) = sum(sum(t_aux))/n_smpl;
        for i = 0:N
            count(k+1,i+1) = count(k+1,i+1) + sum(sum(k_aux == i));
        end
        %I_av(k+1,1) = sum(sum(I_aux))/n_smpl;
        %I_av(k+1,2) = min(I_aux);
        %I_av(k+1,3) = max(I_aux);
    end
end