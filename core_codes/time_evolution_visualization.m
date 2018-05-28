% Time evolution of a system with visualization of the dynamics without
% noise showing the count of nearest neighbors that are ON
close all
clear all
warning off

% lattice parameters
gridsize = 15;
N = gridsize^2;
a0 = 0.5;
Rcell = 0.2*a0;
% circuit parameters
Con = 6;
K = 15;
% initial conditions
p0 = 0.75;
iniON = round(p0*N);
%iniON = 1;

% Initialize parameters
[pos,ex,ey] = init_cellpos_hex(gridsize,gridsize);
dist = dist_mat(pos,gridsize,gridsize,ex,ey);

dist_vec = a0*dist(1,:);
r = dist_vec(dist_vec>0); % exclude self influence
fN = sum(sinh(Rcell)*sum(exp(Rcell-r)./r)); % calculate signaling strength

% generate cell_type (0 case type 1, 1 case type 2)
cell_type = zeros(N,1); % all the same here

% initialize ON cells
cells = zeros(N,1);
cells(randperm(N,iniON)) = 1;
%cells(1:2:N) = 1;

% initialize figure
hin = figure(1);
t = 0;
I = [];
Non = [];
h = [];
cells_hist = {};
update_cell_figure(hin, pos, a0, cells, cell_type, t);
% save vars and update cells
cells_hist{end+1} = cells;
Non(end+1) = sum(cells);
I(end+1) = moranI(cells, a0*dist);
[cells_out, changed, h(end+1)] = update_cells(cells, dist, Con, K, a0, Rcell);
while changed
    pause(1);
    t = t+1;
    %k = waitforbuttonpress;
    update_cell_figure(hin, pos, a0, cells_out, cell_type, t);
    cells_hist{end+1} = cells_out;
    I(end+1) = moranI(cells_out, a0*dist);
    Non(end+1) = sum(cells_out);
    cells = cells_out;
    [cells_out, changed, h(end+1)] = update_cells(cells, dist, Con, K, a0, Rcell);
end
%% Plot trajectory in p, I space

figure();
hold on
plot(Non/N, I, 'r');
plot(Non(1)/N, I(1), 'ro');
plot(Non(end)/N, I(end), 'rx');
xlim([0 1]);
ylim([-0.05 1]);

%% Plot h

%figure();
%plot(0:t, h);
%% Save result
%{
fname_str = sprintf('n%d_neq_%d_a0%d_K_%d_Son_%d_t_%d', ...
    iniON, Non(end), 10*a0, K, Son, t);
i = 1;
fname = fullfile(data_path, 'dynamics_nonoise', ...
    strcat(fname_str,'-v',int2str(i),'.mat'));
while exist(fname, 'file') == 2
    i=i+1;
  fname = fullfile(data_path,'dynamics_nonoise', ...
      strcat(fname_str,'-v',int2str(i),'.mat'));
end
    
save(fname)
%}