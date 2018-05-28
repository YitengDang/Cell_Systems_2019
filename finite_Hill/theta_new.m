% Propose a new spatial order parameter that works for uniform lattices
clear variables
close all

%%
gridsize = 11;
N = gridsize^2;
p = 0.2;
iniON = round(p*N);
a0 = 1.5;

% initiate cells
cells = zeros(N, 1);
cells(randperm(N, iniON)) = 1;
cell_type = zeros(N,1);
%idx = repmat(1:gridsize, 1, round(gridsize)) + kron(0:gridsize:gridsize^2-1, ones(1,gridsize));
idx = repmat(1:gridsize, 1, round(gridsize/2)) + kron(0:2*gridsize:gridsize^2-1, ones(1,gridsize));
%cells(idx) = 1;
%cells = ones(N,1)*0.8786;

% initiate lattice
[pos,ex,ey] = init_cellpos_hex(gridsize,gridsize);
dist = a0*dist_mat(pos,gridsize,gridsize,ex,ey);

% draw cells
h = figure();
update_cell_figure_continuum(h, pos, a0, cells, cell_type, 0);

%% calculate theta
%
cells_pm = 2*cells - 1; % ON cells: 1, OFF cells: -1
cell_mean = mean(cells_pm); % <Xi>

% meshgrid replicates the cells states in x and y directions
[cells_matx, cells_maty] = meshgrid(cells_pm, cells_pm);

% Matrix of cell reading
idx = dist>0;
M = zeros(size(dist));
M(idx) = exp(-dist(idx))./dist(idx);

% Sum of all weights
w_sum = sum(sum(M));
theta = sum(sum(M.*cells_matx.*cells_maty))/w_sum/cell_mean^2;
%}
%% Compare with moranI
%I = moranI(cells, dist);
[I, theta] = moranI(cells, dist);
sprintf('I = %.4f, Theta = %.4f', I, theta)
