function mu_cells = calc_growth_rate(rcell_all, c_growth, K_growth)
    % calculates the growth rate of the cells, according to logistic growth
    A_system = sqrt(3)/2; % area of system
    A_cells = pi*mean(rcell_all.^2); % area covered by cells
    rho = A_cells/A_system; % area density of cells
    mu_cells = c_growth*rho*(1-rho/K_growth);
end