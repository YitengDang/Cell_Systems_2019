function [count, corr_av] = count_eq_parallel_hill_noise_monte_carlo_tmax(...
    dist, Con, K, Rcell, a0, hill, noise, tmax, n_smpl)
% Count the number of equilibrium states found in n_smpl random
% configurations trials, using parallel processing, making the simulation
% faster for large systems.

% This is used to simulate the p_ini vs. p_eq map. It starts with random
% configurations and runs until equilibrium, saving a 2d map, count, that
% for each p_ini and p_eq counts the number of simulations that fall within
% it. t_av counts the average number of steps that it takes to reach
% equilibrium.

% Finishes when tmax is reached, does not check size of fluctuations

% variables
N = size(dist,1); % number of cells
count = zeros(N+1,N+1); % counts of #eq. states
corr_av = zeros(N+1,3);

% simulation parameters
%n_smpl = 1000; % number of samples for each value
%tmax = 1000;
%len_test = 10;

% fix error thresholds
%[err_thr_mean, err_thr_std, ~] = fix_error_thresholds(dist, a0, Con, K, fN, hill, noise);

for k = 0:N % k = <Xi(t=0)> = p_initial
    disp(k)
    corr_aux = zeros(n_smpl,1);
    k_aux = zeros(n_smpl,1);
    %-----------------
    % Test n_smpl random initial configurations and update cells
    % until they are in final configuration
    parfor i = 1:n_smpl
        % initiate cells
        cells = init_random_cells_montecarlo(N, k/N);
        
        % run until tmax reached
        for t=1:tmax
            [cells, ~] = ...
                update_cells_noise_hill(cells, dist, Con, K, a0, Rcell, noise, hill, 10);
        end
        %-----------------
        kout = round(sum(cells)); % round to nearest integer value
        k_aux(i) = kout;
        [~, theta] = moranI(cells, a0*dist);
        corr_aux(i) = theta - (2*mean(cells)-1).^2;
    end
    %-----------------
    % update count
    for i = 0:N
        count(k+1,i+1) = count(k+1,i+1) + sum(sum(round(k_aux) == i));
    end
    % mean, min and max final I
    corr_av(k+1,1) = sum(sum(corr_aux))/n_smpl;
    corr_av(k+1,2) = min(corr_aux);
    corr_av(k+1,3) = max(corr_aux);
end
    
end

