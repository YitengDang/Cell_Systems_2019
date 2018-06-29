% Generates pattern with target I and saves result

% This function generates a pattern with the same fraction of ON cells as
% the vector cells_in and with a spatial order parameter I_target. This is
% done by try and error, turning ON cells that are close to other ON cells.

% Code only works for trying to increase the value of I. For decreasing I
% another code is needed.

close all
clear all
warning off

% lattice parameters
gridsize = 11;
N = gridsize^2;
a0 = 1.5;
Rcell = 0.2*a0;
% circuit parameters
Son = 21;
K = 6;
% initial conditions
p0 = 0.5;
I_target = 0.3;
iniON = round(p0*N);

% Initialize parameters
[pos,ex,ey] = init_cellpos_hex(gridsize,gridsize);
dist = dist_mat(pos,gridsize,gridsize,ex,ey);
dist_vec = a0*dist(1,:);
r = dist_vec(dist_vec>0); % exclude self influence
fN = sum(sinh(Rcell)*sum(exp(Rcell-r)./r)); % calculate signaling strength

% generate cell_type (0 case type 1, 1 case type 2)
cell_type = zeros(N,1); % all the same here

% initialize ON cells
cells_in = zeros(N,1);
cells_in(randperm(N,iniON)) = 1;

% Threshold I 
threshold = 0.001;
%%--------Function code---------------------
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

while abs(I_out-I_target) > threshold && t < 1000
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

disp(I_out);

%% store cells for future use
dir = fullfile(pwd, 'data', 'cells_ini_p_I');
fileid = strrep(sprintf('cells_a0_%.1f_p%.2f_I%.2f', a0, p0, I_out), '.', 'p');

i = 1;
fname = fullfile(dir, strcat(fileid,'-v',int2str(i),'.mat'));
while exist(fname, 'file') == 2
    i=i+1;
  fname = fullfile(dir, strcat(fileid,'-v',int2str(i),'.mat'));
end

save(fname, 'cells_out', 'I_out');
