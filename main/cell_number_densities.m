function [g, gonon, gonoff, dist_vec] = cell_number_densities(cells, dist)
% Calculate the three relevant number densities:
% number density g(r) -> number of cells at a distance r +- dr 
% number of ON cells at a distance r +- dr given cell is ON (gonon(r))
% number of ON cells at a distance r +- dr given cell is OFF (gonoff(r))
% These are calculated based on the distances between cells

N = numel(cells);
Non = sum(cells);
dr = 0.1; % Bin size to calculate the cell number density
dist_vec = dr:dr:18;
num_dist = numel(dist_vec);

% meshgrid replicates the cells states in x and y directions
[~, cells_maty] = meshgrid(cells, cells);
g = zeros(num_dist,1);
gonon = zeros(num_dist,1);
gonoff = zeros(num_dist,1);
for i = 1:num_dist
    idx = dist < dist_vec(i)+dr/2 & dist > dist_vec(i)-dr/2; % get cells in the distance
    g(i) = sum(sum(idx))/2/pi/dist_vec(i)/N/(N-1)/dr; % count cells at the distance
    % only ON cells, first multiplication zeros the line of OFF cells
    gonon(i) = sum((cells_maty.*idx)*cells)/2/pi/dist_vec(i)/N/Non/dr;
    % same as before but the NOT operation gets only the OFF cells
    gonoff(i) = sum(((~cells_maty).*idx)*cells)/2/pi/dist_vec(i)/N/(N-Non)/dr;
end
% 
% r1 = zeros(num_dist,1);
% r1(1:end-1) = (dist_vec(1:end-1) + diff(dist_vec)/2);
% r1(end) = (dist_vec(end) - dist_vec(end-1))/2 + r1(end-1);
% 
% r2 = zeros(num_dist,1);
% r2(2:end) = (dist_vec(2:end) - diff(dist_vec)/2);

% g = g./(r1-r2);
% gonon = gonon./(r1-r2);
% gonoff = gonoff./(r1-r2);

end