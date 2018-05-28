% Time evolution of a system with visualization of the dynamics without
% noise showing the count of nearest neighbors that are ON

close all
clear all
warning off

% Parameters of the system
gridsize = 11;
N = gridsize^2;
a0 = 0.5;
Rcell = 0.2*a0;
iniON = 85;
% parameters
Son = 8;
K = 16;

% Initialize parameters
[pos,ex,ey] = init_cellpos_hex(gridsize,gridsize);
dist = dist_mat(pos,gridsize,gridsize,ex,ey);

% generate cell_type (0 case type 1, 1 case type 2)
cell_type = zeros(N,1); % all the same here

% initialize ON cells
cells = zeros(N,1);
cells(randperm(N,iniON)) = 1;
%cells(1:iniON) = 1;

% initialize figure
hin = figure(1);
t = 0;
update_cell_figure_count_neighbors(hin, pos, dist, a0, cells, cell_type, t);
[cells_out, changed] = update_cells(cells, dist, Son, K, a0, Rcell);

while changed
    t = t+1;
    [~, I(t), ~, ~] = cell_correlation(cells_out, dist);
    update_cell_figure_count_neighbors(hin, pos, dist, a0, cells_out, cell_type, t);
    cells = cells_out;
    [cells_out, changed] = update_cells(cells, dist, Son, K, a0, Rcell);
end

