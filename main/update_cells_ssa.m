function [cells_out, changed, h] = ...
    update_cells_ssa(cells, dist, Son, K, distvals, cvalsN)

% Update cells without noise in a positive feedback loop with infinite hill
% coefficient
% Supplement S.1.1

% Take the SSA result for signaling strength (see hex_lattice_diffusion_ssa.m)

% Account for self-influence
% M: Matrix of cell reading, eq. S10 p.S5
N = numel(cells);
dist2 = round(dist, 6); % crucial: match number of significant digits with number in hex_lattice_diffusion_ssa.m
M = zeros(N, N);
for i=1:N
    for j=1:N
        idx = find(dist2(i,j) == distvals, 1);
        M(i,j) = cvalsN(idx);
    end
end

% Concentration in each cell
C0 = 1 + (Son-1).*cells; % term in brackets of Eq. S9 p.S5

% Reading of each cell
Y = M*C0; % Eq. S9 p.S5

cells_out = Y > K;
changed = ~isequal(cells_out, cells);

% hamiltonian
% h = sum((2*cells-1).*(Y-K)); 
h = sum((cells_out-cells).^2);



        