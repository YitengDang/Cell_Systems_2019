function [count, t_av, I_av, ntmax] = count_eq_parallel_finiteHill_fin_fout(...
    dist, Con, K, Rcell, a0, hill, noise, prec, fp, n_smpl)
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
if numel(fp)~=3
    fprintf('ERROR! Incorrect # fixed points')
    return
end
%%
N = size(dist,1); % number of cells
%n_smpl = 10^3; % number of samples for each value
count = zeros(N+1,N+1); % counts of #eq. states
t_av = zeros(N+1,1);
I_av = zeros(N+1,3);
tmax = 1000;
ntmax = zeros(N+1, 1); % number of simulations reaching tmax
for k = 0:N % k = <Xi(t=0)> = p_initial
    disp(k)
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
        cells = init_random_cells_binary_constraint(N, k, fp);
        [cells, changed] = update_cells_noise_hill(cells, dist, Con, K, a0, Rcell, noise, hill, prec);
        while changed && t_aux(i)<=tmax
            t_aux(i) = t_aux(i) + 1;
            [cells, changed] = update_cells_noise_hill(cells, dist, Con, K, a0, Rcell, noise, hill, prec);
        end
        %--- count tmax---
        if t_aux(i)==tmax
            ntmax_aux(i) = 1;
        end
        %-----------------
        kout = sum(cells>fp(2));
        k_aux(i) = kout;
        [I_aux(i), ~] = moranI(cells, a0*dist);
    end
    %--- count tmax---
    ntmax(k+1) = sum(ntmax_aux);
    %-----------------
    % update count
    for j = 0:N
        count(k+1,j+1) = count(k+1,j+1) + sum(sum(k_aux == j));
    end
    % mean final t
    t_av(k+1) = sum(sum(t_aux))/n_smpl;
    % mean, min and max final I
    I_av(k+1,1) = sum(sum(I_aux))/n_smpl;
    I_av(k+1,2) = min(I_aux);
    I_av(k+1,3) = max(I_aux);
end
    
end

