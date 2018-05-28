function [p_out, I_out, t_av, errors] = count_pI_eq_parallel_2(dist, Con, K, a0, Rcell, p_all, I_all)
% Count the number of equilibrium states found in n_smpl random
% configurations trials, using parallel processing, making the simulation
% faster for large systems.

% This is used to simulate the p_ini vs. p_eq map. It starts with random
% configurations and runs until equilibrium, saving a 2d map, count, that
% for each p_ini and p_eq counts the number of simulations that fall within
% it. t_av counts the average number of steps that it takes to reach
% equilibrium.

N = size(dist,1); % number of cells
n_smpl = 10;
p_out = zeros(numel(p_all), numel(I_all), n_smpl); 
I_out = zeros(numel(p_all), numel(I_all), n_smpl); 
t_av = zeros(numel(p_all), numel(I_all));
dI = (I_all(end)-I_all(1))/(numel(I_all)-1);
errors = 0;

for i=1:numel(p_all)
    p = p_all(i);
    for j=1:numel(I_all)
        I = I_all(j);
        fprintf('p = %.3f, I = %.3f \n', p, I); 
        t_aux = zeros(n_smpl, 1);
        % Test n_smpl random initial configurations and update cells
        % until they are in final configuration
        for k = 1:n_smpl
            fprintf('sample = %d', k);
            % First generate random configuration with required p
            iniON = round(p*N);
            rnd_idx = randperm(N,iniON);
            cells = zeros(N,1);
            cells(rnd_idx) = 1;
            % Then shuffle cells until new I is found
            test = false;
            count = 0;
            cells_in = cells;
            while ~test % generate new lattice with required I
                [cells, ~, ~, test] = generate_I_new(cells_in, I, I+dI, dist, a0);
                count = count+1;
                disp(count);
                if count > 1000
                    disp('Could not generate correct lattice');
                    errors =  errors + 1;
                    break
                end
            end
            % Update cells until equilibrium
            [cells, changed] = update_cells(cells, dist, Con, K, a0, Rcell);
            while changed
                t_aux(k) = t_aux(k) + 1;
                [cells, changed] = update_cells(cells, dist, Con, K, a0, Rcell);
            end
            p_out(i,j,k) = sum(cells)/N;
            [I_out(i,j,k), ~] = moranI(cells, a0*dist); 
        end
        t_av(i,j) = sum(t_aux)/n_smpl;
    end
end

end