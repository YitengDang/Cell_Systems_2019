function [cells_out, changed, mom] = ...
    update_cells_burst(cells, dist, Son, K, a0, Rcell, B)

% Update cells without noise in a positive feedback loop with infinite hill
% coefficient.
% Add a single burst through alteration of the output of a single cell
% B: magnitude of burst, location random site on lattice

% Account for self-influence
% M: Matrix of cell reading, eq. S10 p.S5
idx = dist>0;
M = ones(size(dist)); 
M(idx) = sinh(Rcell)./(a0*dist(idx)).*exp(Rcell-a0*dist(idx));

% Concentration in each cell
if B>0
    %random burst location
    burst_loc = zeros(size(cells));
    burst_loc(randi(size(cells,1))) = 1; 
    C0 = 1 + (Son-1).*cells + B*burst_loc; % term in brackets of Eq. S9 p.S5, with added burst
    %disp(C0);
else
    C0 = 1 + (Son-1).*cells; 
end

% Reading of each cell
Y = M*C0; % Eq. S9 p.S5

mom = sum((2*cells-1).*(Y-K));

cells_out = Y > K;
changed = ~isequal(cells_out, cells);