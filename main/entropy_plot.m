close all
clear all
warning off

% Simulate the entropy map for a set of parameters and plot the data. It
% compares the simulation with the analytical formula

% Parameters of the system
gridsize = 11;
N = gridsize^2;
a0_vec = [1.5]; % all a0 to calculate

% use hexagonal lattice
[dist, pos] = init_dist_hex(gridsize, gridsize);
dist_vec = dist(1,:);

% Runs for all a0 in a0_vec
for idx = 1:numel(a0_vec)
    a0 =a0_vec(idx);
    Rcell = 0.2*a0;
    % --- Entropy as function of Son --- %
    K = 10; % Fixed K
    Son = linspace(1,40,100)'; % Several Son for analytical calculation
    
    % Calculate the analytical result
    omega = zeros(size(Son));
    for i=1:size(Son)
        [omega(i),~] = entropy_eq_sphere(dist_vec, Son(i), K, a0, Rcell);
    end
    
    % Calculate the Monte Carlo simulation for all Son in Son2
    Son2 = linspace(1,40,30)'; % Less points because takes longer
    omega_sim = zeros(size(Son2));
    for i=1:size(Son2)
        [omega_sim(i),~] = sim_entropy(dist, Son2(i), K, a0, Rcell);
        sprintf('Son --> %.0f%%', 100*i/numel(Son2))
    end
    
    % Plot the result
    figure(1)
    plot(Son,omega, 'r-')
    hold on
    plot(Son2,omega_sim, 'ko')
    hold off
    xlabel('S_{ON}')
    ylabel('Entropy')
    % save data
    fname = strcat(datestr(now,'yyyymmdd_hhMM_'), ...
        sprintf('changeSon_N%dK%d_hexagonal.mat', gridsize, K));
    fname = fullfile(pwd, 'data', fname);
    save(fname, 'K', 'Son', 'Son2', 'omega', 'omega_sim' , 'a0', 'Rcell');
    
    % --- Entropy as function of K --- %
    Son = 40; % Fixed Son
    K = linspace(1,80,100)'; % Several K for analytical calculation
    
    % Calculate the analytical result
    omega = zeros(size(K));
    for i=1:size(K)
        [omega(i),~] = entropy_eq_sphere(dist_vec, Son, K(i), a0, Rcell);
    end
    
    % Calculate the Monte Carlo simulation for all K in K2
    K2 = linspace(1,80,30)'; % Less points because takes longer
    omega_sim = zeros(size(K2));
    for i=1:size(K2)
        [omega_sim(i), ~] = sim_entropy(dist, Son, K2(i), a0, Rcell);
        sprintf('K --> %.0f%%', 100*i/numel(K2))
    end
    
    % Plot the result
    figure(2)
    plot(K,omega, 'r-')
    hold on
    plot(K2,omega_sim, 'ko')
    hold off
    xlabel('K')
    ylabel('Entropy')
    % save the data
    fname = strcat(datestr(now,'yyyymmdd_hhMM_'), ...
        sprintf('changeK_N%dSon%d_hexagonal.mat', gridsize, Son));
    fname = fullfile(pwd, 'data', fname);
    save(fname, 'K', 'Son', 'K2', 'omega', 'omega_sim', 'a0', 'Rcell');
end
