%% Try to find lattice configurations with desired energies by brute force sampling
h_target = 12:2:22;
margin = 1;
target_cells = zeros(length(h_target), L^2);

for i=1:length(h_target)
    found = 0;
    while ~found
        cells = randi(2,L^2,1)-1; %random lattice
        h_rand = hamiltonian(cells, dist, Son, K, a0, Rcell); % energy of lattice
        if abs(h_rand-h_target(i)) < margin
            target_cells(i,:) = cells;
            disp(h_target(i));
            disp(cells);
            found = 1;
        end
    end
end

%}