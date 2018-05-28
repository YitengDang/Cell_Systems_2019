function [count, t_av] = count_eq_parallelbatch(dist, Son, K, a0, Rcell)
% Count the number of equilibrium states found in n_smpl ran,dom
% configurations trials

N = size(dist,1); % number of cells
n_smpl = 250;
count = zeros(N+1,N+1);
t_av = zeros(N+1,1);
for k = 0 : N
    if k == 0
        cells = zeros(N,1); % initialize all cells OFF
        [cells, changed] = update_cells(cells, dist, Son, K, a0, Rcell);
        while changed
            t_av(k+1) = t_av(k+1) + 1;
            [cells, changed] = update_cells(cells, dist, Son, K, a0, Rcell);
        end
        kout = sum(cells);
        count(k+1,kout+1) = count(k+1,kout+1) + 1;
    else
        if k == N
            cells = ones(N,1); % initialize all cells ON
            [cells, changed] = update_cells(cells, dist, Son, K, a0, Rcell);
            while changed
                t_av(k+1) = t_av(k+1) + 1;
                [cells, changed] = update_cells(cells, dist, Son, K, a0, Rcell);
            end
            kout = sum(cells);
            count(k+1,kout+1) = count(k+1,kout+1) + 1;
        else
            % Test n_smpl different sampling and see if it is eq.
            for i = 1:n_smpl
                rnd_idx = randperm(N,k);
                cells = zeros(N,1);
                cells(rnd_idx) = 1;
                [cells, changed] = update_cells(cells, dist, Son, K, a0, Rcell);
                while changed
                    t_av(k+1) = t_av(k+1) + 1;
                    [cells, changed] = update_cells(cells, dist, Son, K, a0, Rcell);
                end
                kout = sum(cells);
                count(k+1,kout+1) = count(k+1,kout+1) + 1;
            end
            t_av(k+1) = t_av(k+1)/n_smpl;
        end
    end
end

end

