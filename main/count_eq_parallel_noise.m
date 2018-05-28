function [count, t_av, I_av] = count_eq_parallel_noise(dist, Son, K, a0, Rcell, noise)
% Count the number of equilibrium states found in n_smpl random
% configurations trials, using parallel processing, making the simulation
% faster for large systems. This function assumes that there is noise in
% the system and measures the stop condition by doing a moving std of
% the number of ON cells and the spatial order parameter, stopping if the
% std of this window is below the threshold specified.

% This is used to simulate the p_ini vs. p_eq map. It starts with random
% configurations and run until equilibrium, saving a 2d map, count, that
% for each p_ini and p_eq counts the number of simulations that fall within
% it. t_av counts the average number of steps that it takes to reach
% equilibrium. I_av is the average spatial order parameter among the
% simulations

N = size(dist,1); % number of cells
n_smpl = 1000; % number of samples to generate the statistics

% These variables are to stop the running
err_thr_n = 1; % variation allowed for number of ON cells
err_thr_I = 0.02; % variation allowed for I
len_test = 10; % window (in steps) of the moving std
t_max = 1000; % if simulation exceeds t_max, do not register result

count = zeros(N+1,N+1);
t_av = zeros(N+1,1);
I_av = zeros(N+1,3);
for k = 0 : N
    disp(k)
    t_aux = zeros(n_smpl,1);
    k_aux = zeros(n_smpl,1);
    I_aux = zeros(n_smpl,1);
    parfor i = 1:n_smpl
        % initialize cells
        if k == 0
            cells = zeros(N,1); % initialize all cells OFF
        else
            if k == N
                cells = ones(N,1); % initialize all cells ON
            else
                % Test n_smpl different sampling and see if it is eq.
                rnd_idx = randperm(N,k);
                cells = zeros(N,1);
                cells(rnd_idx) = 1;
            end
        end
        cont = true;
        I = zeros(len_test,1);
        fracON = zeros(len_test,1);
        t = 0;
        while cont
            % First run to populate the moving average test
            if t < len_test
                t = t+1;
                ind_aux = len_test - t+1; % populate in reverse order (older last)
                [cells, ~, ~, I(ind_aux)] = ...
                    update_cells_noise(cells, dist, Son, K, a0, Rcell, noise);
                fracON(ind_aux) = sum(cells);
            elseif t > t_max
                disp('DNF; t_max reached');
                break
            else
                % Check if the std falls within the threshold
                if std(fracON) < err_thr_n && std(I) < err_thr_I
                    cont = false;
                    t_aux(i) = t;
                    k_aux(i) = round(mean(fracON));
                    I_aux(i) = mean(I);
                else
                    t = t+1;
                    I = circshift(I,1); % rotate the list to replace the last 
                    fracON = circshift(fracON, 1);
                    [cells, ~, ~, I(1)] = ...
                        update_cells_noise(cells, dist, Son, K, a0, Rcell, noise);
                    fracON(1) = sum(cells);
                end
            end
        end
    end
    t_av(k+1) = sum(sum(t_aux))/n_smpl;
    for n = 0:N
        count(k+1,n+1) = count(k+1,n+1) + sum(sum(k_aux == n));
    end
    I_av(k+1,1) = sum(sum(I_aux))/n_smpl;
    I_av(k+1,2) = min(I_aux);
    I_av(k+1,3) = max(I_aux);
end

end
