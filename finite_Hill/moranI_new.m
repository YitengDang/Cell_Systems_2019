function [I2, theta] = moranI_new(cells, dist)
% This function calculates a new correlation function that is more
% appropriate for uniform lattices
% The weights are assumed to be the interaction f_ij
% Enter dist in correct units (of a0)

cells_pm = 2*cells - 1; % ON cells: 1, OFF cells: -1
% meshgrid replicates the cells states in x and y directions
[cells_matx, cells_maty] = meshgrid(cells_pm, cells_pm);

% For the Moran I we use the factor of the diffusion equation as weight
% Account for self-influence
idx = dist>0;
M = zeros(size(dist));

% Matrix of cell reading
M(idx) = exp(-dist(idx))./dist(idx);

% Sum of all weights
w_sum = sum(sum(M));

% Calculate Theta and I2
theta = sum(sum(M.*cells_matx.*cells_maty)); % not normalized by w_sum
I2 = theta - w_sum*mean(cells_pm)^2;