function [cells_out, changed, dq, dw, H22] = ...
    update_cells_tune_B_field_varya0(cells, dist, a0, Rcellf, Son, K, a0_old)

% Update cells for a system with changing B field
% System without noise in a positive feedback loop with infinite hill
% coefficient
% Supplement S.1.1

% Account for self-influence
% M: Matrix of cell readings, eq. S10 p.S5
idx = dist>0;
M1 = ones(size(dist));
M2 = ones(size(dist));
M1(idx) = sinh(a0_old*Rcellf)./(a0_old*dist(idx)).*exp(a0_old*Rcellf-a0_old*dist(idx));
M2(idx) = sinh(a0*Rcellf)./(a0*dist(idx)).*exp(a0*Rcellf-a0*dist(idx));

% Concentration in each cell
C0 = 1 + (Son-1).*cells; % term in brackets of Eq. S9 p.S5

% Reading of each cell
% Eq. S9 p.S5
Y11 = M1*C0;
Y12 = M2*C0; 

% update cells
cells_out = Y12 > K;
changed = ~isequal(cells_out, cells);

% new readings
C0_22 = 1 + (Son-1).*cells_out;
Y22 = M2*C0_22;

% energies
H11 = -sum((2*cells-1).*(Y11-K));
H12 = -sum((2*cells-1).*(Y12-K)); 
H22 = -sum((2*cells_out-1).*(Y22-K)); 

% energy changes
N = numel(cells);
dw = (H12 - H11)/N;
dq = (H22 - H12)/N;

%-----Calculate work------
%{
% Old method, did not work (contains error)
% magnetic work
dB = ((Son+1)/2*(1+fN)-K)-((Son_old+1)/2*(1+fN_old)-K_old);
dwB = -dB*(2*sum(cells)/numel(cells)-1); % dW = - dB (2p-1)

% remaining work
M_old = ones(size(dist)); % old matrix of f(r_{ij})
M_old(idx) = sinh(Rcell)./(a0_old*dist(idx)).*exp(Rcell-a0_old*dist(idx)); 

cells_pm = 2*cells - 1; % ON cells: 1, OFF cells: -1
[cellsX, cellsY] = meshgrid(cells_pm, cells_pm);
theta = sum(sum(M.*cellsX.*cellsY));
theta_old = sum(sum(M_old.*cellsX.*cellsY));
dwR = - (Son-1)/2 * theta + (Son_old-1)/2 * theta_old - (Son - Son_old)/2;  
%}