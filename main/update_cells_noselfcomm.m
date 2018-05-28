function [cells_out, changed, mom] = ...
    update_cells_noselfcomm(cells, dist, Son, K, a0, Rcell)

% Update cells without noise in a positive feedback loop with infinite hill
% coefficient

% Account for self-influence
idx = dist>0;
M = zeros(size(dist));

% Matrix of cell reading
M(idx) = sinh(Rcell)./(a0*dist(idx)).*exp(Rcell-a0*dist(idx));

% Concentration in each cell
C0 = 1 + (Son-1).*cells;

% Reading of each cell
Y = M*C0;

mom = sum((2*cells-1).*(Y-K));

cells_out = Y > K;
changed = ~isequal(cells_out, cells);



        