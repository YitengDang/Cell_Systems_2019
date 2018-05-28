% Calculate the correlations for a predetermined configuration. It
% calculates the number density and the correlation <Xi Xj> where Xi is the
% cell state of cell i.

close all
clear variables
warning off

% Parameters of the system
gridsize = 21;
N = gridsize^2;
a0 = 0.5;
Rcell = 0.2*a0;
save_fig = 1;

% initialize hexagonal grid
[pos,ex,ey] = init_cellpos_hex(gridsize,gridsize);
dist = dist_mat(pos,gridsize,gridsize,ex,ey);

% Define the configuration to be tested
cells = zeros(N,1);
Non = 50;
%cells(randperm(N,Non)) = 1;
cells(1:Non) = 1;

% generate cell_type (0 case type 1, 1 case type 2)
cell_type = zeros(N,1); % all the same here

% Calculate the cell correlation in space and the spatial order parameter
[cc, I, ~, r] = cell_correlation(cells, dist);

% Plot the cell configuration
hin = figure(1);
t = 0;
update_cell_figure(hin, pos, a0, cells, cell_type, t);

% Plot the correlation in distance
figure(2)
plot(r(2:end), cc(2:end))
ylabel('<X_i X_j>')
xlabel('r')

% Plot the cell number densities for ON cells and the average cell number
% density
[g, gonon, gonoff, r] = cell_number_densities(cells, dist);
figure(3)
gon = ((N-Non)*gonoff+ Non*gonon)/N;
plot(r, gon, r, g)
legend({'ON cells';'All cells'})
ylabel('Cell number density')
xlabel('r')


