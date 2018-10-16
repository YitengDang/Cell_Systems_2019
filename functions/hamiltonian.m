function h=hamiltonian(cells, dist, Son, K, a0, Rcell)
    % Calculates the Hamiltonian of a given lattice of cells
    %disp(lattice);
    
    idx = dist>0;
    M = ones(size(dist)); 
    M(idx) = sinh(Rcell)./(a0*dist(idx)).*exp(Rcell-a0*dist(idx));

    % Concentration in each cell
    C0 = 1 + (Son-1).*cells; % term in brackets of Eq. S9 p.S5

    % Reading of each cell
    Y = M*C0; % Eq. S9 p.S5
    
    % Hamiltonian
    h = sum(sum(-cells.*(Y-K)));
end