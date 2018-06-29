function [cells_out, changed, dq, dw, H22] = ...
    update_cells_tune_B_field(cells, dist, a0, Rcell, Son, K, Son_old, K_old)

% Update cells for a system with changing B field
% System without noise in a positive feedback loop with infinite hill
% coefficient
% Supplement S.1.1
%
% Account for self-influence
% M: Matrix of cell readings, eq. S10 p.S5
idx = dist>0;
M = ones(size(dist));
M(idx) = sinh(Rcell)./(a0*dist(idx)).*exp(Rcell-a0*dist(idx));

% Concentration in each cell
C0_12 = 1 + (Son-1).*cells; % term in brackets of Eq. S9 p.S5
C0_11 = 1 + (Son_old-1).*cells; % term in brackets of Eq. S9 p.S5

% Reading of each cell
Y12 = M*C0_12; % Eq. S9 p.S5
Y11 = M*C0_11;

% update cells
cells_out = Y12 > K;
changed = ~isequal(cells_out, cells);

% new readings
C0_22 = 1 + (Son-1).*cells_out;
Y22 = M*C0_22;

% energies
H11 = -sum((2*cells-1).*(Y11-K_old));
H12 = -sum((2*cells-1).*(Y12-K)); 
H22 = -sum((2*cells_out-1).*(Y22-K)); 

% energy changes
N = numel(cells);
dw = (H12 - H11)/N;
dq = (H22 - H12)/N;

%-----Calculate work------
% to separate file
%{
% magnetic work
dB = ((Son+1)/2*(1+fN)-K)-((Son_old+1)/2*(1+fN)-K_old);
dwB = -dB*(2*sum(cells)/numel(cells)-1); % dw = - dB (2p-1)

% remaining work
cells_pm = 2*cells - 1; % ON cells: 1, OFF cells: -1
[cellsX, cellsY] = meshgrid(cells_pm, cells_pm);
theta = sum(sum(M.*cellsX.*cellsY))/numel(cells); % theta = Theta/N
dwR = - (Son-Son_old)/2 * (theta+1);  
%}