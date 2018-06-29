function [cells_out, I_out, t] = generate_pattern(cells_in, I_target, dist)
% This function generates a pattern with the same fraction of ON cells as
% the vector cells_in and with a spatial order parameter I_target. This is
% done by try and error, turning ON cells that are close to other ON cells.

% Get the number of cells
N = numel(cells_in);
% vector of indices
idx_vec = transpose(1:N);

% First neighbor distance
eps = 1e-5;
dist_vec = get_unique_distances(dist, eps);
dist1 = dist_vec(2);
first_nei = 1*(dist < dist1+eps & dist > dist1-eps);

% Calculate the initial spatial order
I_out = moranI(cells_in, dist);
cells_out = cells_in;
p = sum(cells_in)/N;
t = 0;
while I_out < I_target && t < 1000
    t = t+1;
    % number of neighbors that are ON
    nei_ON = first_nei*cells_out;
    
    % Get the indexes of the ON cells with at least one OFF neighbor
    idx_ON = idx_vec((cells_out == 1) & (nei_ON < 6));
    
    if numel(idx_ON) == 0
        cells_out = cells_in;
        I_out = moranI(cells_out, dist);
        break
    end
    % Randomly pick an ON cell
    idx = datasample(idx_ON, 1);

    % Find the neighbors
    neighbors = idx_vec(first_nei(idx,:) == 1);

    % Turn neighbors ON
    cells_out(neighbors) = 1;
    
    % Turn ON cell OFF. Before doing so remove the set already used from
    % the list.
    idx_ON = idx_vec((cells_out == 1) & (nei_ON < 6*p+1));
    idx_ON(ismember(idx_ON, [neighbors; idx])) = [];
    if numel(idx_ON) == 0
        cells_out = cells_in;
        I_out = moranI(cells_out, dist);
        break
    end
    cells_out(datasample(idx_ON, 6-nei_ON(idx),'Replace',false)) = 0;

    I_out = moranI(cells_out, dist);  
end
        