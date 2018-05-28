function [count, t_av, corr_av, ntmax] = count_eq_parallel_hill_noise_monte_carlo(...
    dist, Con, K, Rcell, a0, fN, hill, noise)
% Count the number of equilibrium states found in n_smpl random
% configurations trials, using parallel processing, making the simulation
% faster for large systems.

% This is used to simulate the p_ini vs. p_eq map. It starts with random
% configurations and runs until equilibrium, saving a 2d map, count, that
% for each p_ini and p_eq counts the number of simulations that fall within
% it. t_av counts the average number of steps that it takes to reach
% equilibrium.

% variables
N = size(dist,1); % number of cells
count = zeros(N+1,N+1); % counts of #eq. states
t_av = zeros(N+1,1);
corr_av = zeros(N+1,3);
ntmax = zeros(N+1, 1); % number of simulations reaching tmax

% simulation parameters
n_smpl = 1000; % number of samples for each value
tmax = 2000;
len_test = 10;

% fix error thresholds
[err_thr_mean, err_thr_std, ~] = fix_error_thresholds(dist, a0, Con, K, fN, hill, noise);

for k = 0:N % k = <Xi(t=0)> = p_initial
    disp(k)
    ntmax_aux = zeros(n_smpl,1);
    corr_aux = zeros(n_smpl,1);
    k_aux = zeros(n_smpl,1);
    t_aux = zeros(n_smpl,1);
    %-----------------
    % Test n_smpl random initial configurations and update cells
    % until they are in final configuration
    parfor i = 1:n_smpl
        % initiate cells
        cells = init_random_cells_montecarlo(N, k/N);
        
        % first populate the test
        xmean_last = zeros(len_test, 1);
        xstd_last = xmean_last;
        for t=1:len_test
            ind_aux = len_test - t+1; % populate in reverse order (older last)
            [cells, ~] = ...
                update_cells_noise_hill(cells, dist, Con, K, a0, Rcell, noise, hill, 10);
            % update test variables
            xmean_last(ind_aux) = mean(cells);
            xstd_last(ind_aux) = std(cells);
        end
        
        % run until errors below thresholds
        cont = true;
        while cont
            if std(xmean_last) < err_thr_mean && std(xstd_last) < err_thr_std
                cont = false;
            elseif t > tmax
                disp('DNF; t_max reached');
                cont = false;
                ntmax_aux(i) = 1;
            else
                t = t+1;
                [cells, ~] = ...
                    update_cells_noise_hill(cells, dist, Con, K, a0, Rcell, noise, hill, 10);
                % rotate the list to replace the last 
                xmean_last = circshift(xmean_last, 1);
                xstd_last = circshift(xstd_last, 1);
                % update test variables
                xmean_last(1) = mean(cells);
                xstd_last(1) = std(cells);
            end
        end
        %-----------------
        t_aux(i) = t;
        kout = round(sum(cells)); % round to nearest integer value
        k_aux(i) = kout;
        [~, theta] = moranI(cells, a0*dist);
        corr_aux(i) = theta - (2*mean(cells)-1).^2;
    end
    %--- count ntmax---
    ntmax(k+1) = sum(ntmax_aux);
    %-----------------
    % update count
    for i = 0:N
        count(k+1,i+1) = count(k+1,i+1) + sum(sum(round(k_aux) == i));
    end
    % mean final t
    t_av(k+1) = sum(sum(t_aux))/n_smpl;
    % mean, min and max final I
    corr_av(k+1,1) = sum(sum(corr_aux))/n_smpl;
    corr_av(k+1,2) = min(corr_aux);
    corr_av(k+1,3) = max(corr_aux);
end
    
end

