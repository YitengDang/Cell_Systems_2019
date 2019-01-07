function [rcell_all_new, delta_A_cells] = calc_cell_sizes_logistic_growth(rcell_all, c_growth, K_growth)
    % calculates the new sizes of the cells, according to logistic growth
    % of each individual cell
    if c_growth == 0
        rcell_all_new = rcell_all;
        delta_A_cells = 0;
        return
    end
    A_system = sqrt(3)/2; % area of system
    A_cells = pi*rcell_all.^2; % area covered by cells
    rho = A_cells./A_system; % area density of cells
    rcell_all_new = c_growth.*rho.*(1+1/c_growth-rho./K_growth);
    
    % Cells that attain a negative radius are set to a very small size
    rcell_all_new(rcell_all_new<0) = 0.01;
    
    % Determine change in area of all cells
    delta_A_cells = pi*(rcell_all_new.^2 - rcell_all.^2);
end