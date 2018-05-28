function [dW] = tune_B_calculate_work(cells, dist, a0, Rcell, Son, K, Son_old, K_old)
% M: Matrix of cell readings, eq. S10 p.S5
idx = dist>0;
M = ones(size(dist));
M(idx) = sinh(Rcell)./(a0*dist(idx)).*exp(Rcell-a0*dist(idx));

% Concentration in each cell
C0 = 1 + (Son-1).*cells; % term in brackets of Eq. S9 p.S5
C0_old = 1 + (Son_old-1).*cells; % term in brackets of Eq. S9 p.S5

% Reading of each cell
Y = M*C0; % Eq. S9 p.S5
Y_old = M*C0_old;

% energies
H = -sum((2*cells-1).*(Y-K)); %-energy of the new system 
H_old = -sum((2*cells-1).*(Y_old-K_old)); %-energy of the old system
dW = (H - H_old)/numel(cells);

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
end