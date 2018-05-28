close all
clear all
warning off

% Simulate the entropy map for a set of parameters and plot the data

% Parameters of the system
gridsize = 11;
N = gridsize^2;
a0_vec = [0.5; 1; 1.5];

% use hexagonal lattice
[pos,ex,ey] = init_cellpos_hex(gridsize,gridsize);
dist = dist_mat(pos,gridsize,gridsize,ex,ey);
dist_vec = dist(1,:);

for idx = 1:numel(a0_vec)
    a0 =a0_vec(idx);
    Rcell = 0.2*a0;
    % --- Entropy as function of Son --- %
    K = 10;
    Son = linspace(1,40,100)';
    
    omega = zeros(size(Son));
    for i=1:size(Son)
        [omega(i),~] = entropy_eq_sphere(dist_vec, Son(i), K, a0, Rcell);
    end
    
    Son2 = linspace(1,40,30)';
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
    Son = 40;
    K = linspace(1,80,100)';
    
    omega = zeros(size(K));
    for i=1:size(K)
        [omega(i),~] = entropy_eq_sphere(dist_vec, Son, K(i), a0, Rcell);
    end
    
    K2 = linspace(1,80,30)';
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
    fname = strcat(datestr(now,'yyyymmdd_hhMM_'), ...
        sprintf('changeK_N%dSon%d_hexagonal.mat', gridsize, Son));
    fname = fullfile(pwd, 'data', fname);
    save(fname, 'K', 'Son', 'K2', 'omega', 'omega_sim', 'a0', 'Rcell');
end
