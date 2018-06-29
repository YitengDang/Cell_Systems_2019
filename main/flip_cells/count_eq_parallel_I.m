function [count, t_av] = count_eq_parallel_I(dist, Son, K, a0, Rcell, I)
% Count the number of equilibrium states found in n_smpl random
% configurations trials, using parallel processing, making the simulation
% faster for large systems. In this function, instead of having a random
% configuration, it tries to generate a random pattern that has an I close
% to the specified value I that is input to the function.

% This is used to simulate the p_ini vs. p_eq map. It starts with random
% configurations and run until equilibrium, saving a 2d map, count, that
% for each p_ini and p_eq counts the number of simulations that fall within
% it. t_av counts the average number of steps that it takes to reach
% equilibrium.

% Used in: initial_vs_eq_countON_I_batch_save.m and
% initial_vs_eq_countON_I_plotdata.m files

N = size(dist,1); % number of cells
% We make the graph from p = 0.2 to p = 0.8
N1 = round(0.2*N);
N2 = round(0.8*N);
n_smpl = 500;
count = zeros(N2-N1+1,N+1);
t_av = zeros(N2-N1+1,1);
for k = N1 : N2
    disp(k)
    
    % Test n_smpl different sampling and see if it is eq.
    t_aux = zeros(n_smpl,1);
    k_aux = zeros(n_smpl,1);
    worked = zeros(n_smpl,1);
    parfor i = 1:n_smpl
        cells = zeros(N,1);
        cells(randperm(N, k)) = 1;
        % Generate a pattern with the desired I. Only counts the trials
        % that worked.
        [cells, I_out] = generate_pattern(cells, I, a0*dist);
        if I_out >= I
            worked(i) = 1;
            [cells, changed] = update_cells(cells, dist, Son, K, a0, Rcell);
            while changed
                t_aux(i) = t_aux(i) + 1;
                [cells, changed] = update_cells(cells, dist, Son, K, a0, Rcell);
            end
            kout = sum(cells);
            k_aux(i) = kout;
        end
    end
    t_av(k-N1+1) = sum(sum(t_aux))/sum(worked);
    for i = 0:N
        count(k-N1+1,i+1) = count(k-N1+1,i+1) + sum(sum(k_aux == i));
    end
end

end
