clear variables
close all

% Parameters

I_target = 0.2;
grid = 15;
N = grid^2;
p = 0.8;
a0 = 0.5;

[dist, pos] = init_dist_hex(grid, grid);

cells = zeros(N,1);
cells(randperm(N,round(p*N))) = 1;

hini = figure(1);
update_cell_figure_withI(hini, pos, dist, a0, cells, zeros(N, 1), 0)

[cells_out, I_out] = generate_pattern(cells, I_target, a0*dist);
hend = figure(2);
update_cell_figure_withI(hend, pos, dist, a0, cells_out, zeros(N, 1), 0)