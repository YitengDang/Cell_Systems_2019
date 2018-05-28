close all
clear all
warning off

% Test the correlations for different configurations

% Parameters of the system
gridsize = 21;
N = gridsize^2;
a0 = 0.5;
Rcell = 0.2*a0;
save_fig = 1;

% use hexagonal lattice
[pos,ex,ey] = init_cellpos_hex(gridsize,gridsize);
dist = dist_mat(pos,gridsize,gridsize,ex,ey);
cells = zeros(N,1);
Non = 50;
%cells(randperm(N,Non)) = 1;
cells(1:Non) = 1;

% generate cell_type (0 case type 1, 1 case type 2)
cell_type = zeros(N,1); % all the same here

[cc, I, r] = cell_correlation(cells, dist);

% initialize figure
hin = figure(1);
t = 0;
update_cell_figure(hin, pos, a0, cells, cell_type, t);

figure(2)
plot(r(2:end), cc(2:end))

[g, gonon, gonoff, r] = cell_number_densities(cells, dist);
figure(3)
gon = ((N-Non)*gonoff+ Non*gonon)/N;
plot(r, gon, r, g)


