function [cc, I, nn, dist_vec] = cell_correlation(cells, dist, dist_vec)
% Calculate the correlation of cells, the spatial correlation
% measured by the Moran Index and the average of ON nearest neighbors

eps = 1e-4;

if nargin < 3
    dist_vec = get_unique_distances(dist, eps);
end

cells_pm = 2*cells - 1; % ON cells: 1, OFF cells: -1
cell_mean = mean(cells_pm); % <Xi>

% First neighbor distance
dist1 = dist_vec(2);
idx = 1*(dist < dist1+eps & dist > dist1-eps);
nn = mean(idx*cells);

% multiply the state of every cell to those at a certain distance
% meshgrid replicates the cells states in x and y directions
[cells_matx, cells_maty] = meshgrid(cells_pm, cells_pm);
cc = zeros(numel(dist_vec),1);
for i = 1:numel(dist_vec)
    idx = dist < dist_vec(i)+eps & dist > dist_vec(i)-eps; % get cells in the distance
    cells_mat_mult = cells_maty.*idx;
    tmp = (cells_mat_mult*cells_pm)./sum(idx,2);
    cc(i) = mean(tmp) - cell_mean^2; % corr = <Xi Xj> - <Xi><Xj>
end

% For the Moran I we use the factor of the diffusion equation as weight
% Account for self-influence
idx = dist>0;
M = zeros(size(dist));

% Matrix of cell reading
M(idx) = exp(-dist(idx))./dist(idx);

% Sum of all weights
w_sum = sum(sum(M));

% cells state variance (using N as normalization)
cells_var = var(cells_pm,1);

% Multiply every weight by the cell state in a crossed fashion
tmp = M.*(cells_matx-cell_mean).*(cells_maty-cell_mean);
I = sum(sum(tmp))/w_sum/cells_var;

end
            
        
        
