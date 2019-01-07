function [cells_out, changed] = ...
    update_cells_survival(cells, dist, S, K, Rcell)
% Works for both closed and periodic bcs. Boundary effects completely go
% into dist.

% M: Matrix of cell reading
idx = dist>0;
M = ones(size(dist)); 
M(idx) = sinh(Rcell)./(dist(idx)).*exp(Rcell-dist(idx));

% Secretion rate of each cell
C0 = S*cells; 

% Reading of each cell
Y = M*C0; 

cells_out = Y > K;
changed = ~isequal(cells_out, cells);

end