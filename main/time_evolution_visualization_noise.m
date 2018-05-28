% Time evolution of a system with visualization of the dynamics with noise

close all
clear all
warning off

% Parameters of the system
gridsize = 21;
N = gridsize^2;
a0 = 1;
Rcell = 0.2*a0;
iniON = round(0.5*N);
eq_before = 0;
load_last = 0;
t_max = 200;
saveq = 0; %0: don't save, 1: save

% parameters
Son = 12;
K = 13;
noise = 0.15*K;

% Initialize parameters
[pos,ex,ey] = init_cellpos_hex(gridsize,gridsize);
dist = dist_mat(pos,gridsize,gridsize,ex,ey);

dist_vec = dist(1,:);
r = a0*dist_vec(dist_vec>0); % exclude self influence
fN = sum(sinh(Rcell)*sum(exp(Rcell-r)./r));

% generate cell_type (0 case type 1, 1 case type 2)
cell_type = zeros(N,1); % all the same here

% initialize ON cells
if load_last > 0
    load('last_cells.mat')
else
    cells = zeros(N,1);
    cells(randperm(N,iniON)) = 1;
    %cells(1:iniON) = 1;
end

if eq_before > 0
    changed = true;
    while changed
        [cells_out, changed, ~] = update_cells(cells, dist, Son, K, a0, Rcell);
        cells = cells_out;
    end
    save('last_cells.mat', 'cells')
end
I = zeros(t_max+1,1);
Non = zeros(t_max+1,1);
mom = zeros(t_max+1,1);
cells_hist = {};

% initialize figure
hin = figure(1);
t = 0;
update_cell_figure(hin, pos, a0, cells, cell_type, t);
cells_hist{end+1} = cells;
F(1) = getframe(gcf);
Non(t+1) = sum(cells);
I(t+1) = moranI(cells, a0*dist);
[cells_out, changed, mom(t+1)] = update_cells_noise(cells, dist, Son, K, a0, Rcell, noise);

for t = 1:t_max
    %pause(1);
    update_cell_figure(hin, pos, a0, cells_out, cell_type, t);
    cells_hist{end+1} = cells_out;
    I(t+1) = moranI(cells_out, a0*dist);
    Non(t+1) = sum(cells_out);
    F(t+1) = getframe(gcf);
    cells = cells_out;
    [cells_out, changed, mom(t+1)] = update_cells_noise(cells, dist, Son, K, a0, Rcell, noise);
end

fname_str = sprintf('N%d_n%d_neq_%d_a0%d_K_%d_Con_%d_noise_%d', ...
    N, iniON, Non(end), round(10*a0), round(K), round(Son), round(10*noise));
i = 1;
fname = fullfile(pwd, 'data', 'dynamics_noise', ...
    strcat(fname_str,'-v',int2str(i),'.mat'));
while exist(fname, 'file') == 2
    i=i+1;
    fname = fullfile(pwd,'data','dynamics_noise', ...
        strcat(fname_str,'-v',int2str(i),'.mat'));
end

if saveq
    save(fname)
end