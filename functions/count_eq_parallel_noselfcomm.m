function [count, t_av, I_av] = count_eq_parallel_noselfcomm(dist, Son, K, a0, Rcell)
% Count the number of equilibrium states found in n_smpl random
% configurations trials, using parallel processing, making the simulation
% faster for large systems.

% This is used to simulate the p_ini vs. p_eq map. It starts with random
% configurations and run until equilibrium, saving a 2d map, count, that
% for each p_ini and p_eq counts the number of simulations that fall within
% it. t_av counts the average number of steps that it takes to reach
% equilibrium.

N = size(dist,1); % number of cells
n_smpl = 1000;
count = zeros(N+1,N+1);
t_av = zeros(N+1,1);
I_av = zeros(N+1,3);
for k = round(0.4*N) : round(0.6*N)
    disp(k)
    if k == 0
        cells = zeros(N,1); % initialize all cells OFF
        [cells, changed] = update_cells_noselfcomm(cells, dist, Son, K, a0, Rcell);
        while changed
            t_av(k+1) = t_av(k+1) + 1;
            [cells, changed] = update_cells_noselfcomm(cells, dist, Son, K, a0, Rcell);
        end
        kout = sum(cells);
        count(k+1,kout+1) = count(k+1,kout+1) + 1;
    else
        if k == N
            cells = ones(N,1); % initialize all cells ON
            [cells, changed] = update_cells_noselfcomm(cells, dist, Son, K, a0, Rcell);
            while changed
                t_av(k+1) = t_av(k+1) + 1;
                [cells, changed] = update_cells_noselfcomm(cells, dist, Son, K, a0, Rcell);
            end
            kout = sum(cells);
            count(k+1,kout+1) = count(k+1,kout+1) + 1;
        else
            % Test n_smpl different sampling and see if it is eq.
            t_aux = zeros(n_smpl,1);
            I_aux = zeros(n_smpl,1);
            k_aux = zeros(n_smpl,1);
            for i = 1:n_smpl
                rnd_idx = randperm(N,k);
                cells = zeros(N,1);
                cells(rnd_idx) = 1;
                cells_ini = cells;
                [cells, changed] = update_cells_noselfcomm(cells, dist, Son, K, a0, Rcell);
                while changed
                    t_aux(i) = t_aux(i) + 1;
                    [cells, changed] = update_cells_noselfcomm(cells, dist, Son, K, a0, Rcell);
                    if t_aux(i) > 1000
                        save('tmp_data.mat', 'cells_ini');
                        changed = false;
                    end
                end
                kout = sum(cells);
                k_aux(i) = kout;
                [I_aux(i), ~] = moranI(cells, a0*dist);
            end
            t_av(k+1) = sum(sum(t_aux))/n_smpl;
            for i = 0:N
                count(k+1,i+1) = count(k+1,i+1) + sum(sum(k_aux == i));
            end
            I_av(k+1,1) = sum(sum(I_aux))/n_smpl;
            I_av(k+1,2) = min(I_aux);
            I_av(k+1,3) = max(I_aux);
        end
    end
end

end
