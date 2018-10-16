function [cells_out, changed, h, hlist] = ...
    update_cells(cells, dist, Con, K, a0, Rcell)
%%
% Update cells without noise in a positive feedback loop with infinite hill
% coefficient
% Supplement S.1.1

% Account for self-influence
% M: Matrix of cell reading, eq. S10 p.S5
idx = dist>0;
M = ones(size(dist)); 
M(idx) = sinh(Rcell)./(a0*dist(idx)).*exp(Rcell-a0*dist(idx));

% Concentration in each cell
C0 = 1 + (Con-1).*cells; % term in brackets of Eq. S9 p.S5

% Reading of each cell
Y = M*C0; % Eq. S9 p.S5

cells_out = Y > K;
changed = ~isequal(cells_out, cells);

% hamiltonian
hlist = -(2*cells-1).*(Y-K);
h = sum(hlist);

        