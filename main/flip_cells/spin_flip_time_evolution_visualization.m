% Generate trajectories of similar initial conditions by switching single
% cells 
% Visualizes effect of flipping one or more cells of an equilibrium
% configuration

% Time evolution of a system without noise and without visualization
close all
clear all
warning off

% Parameters of the system
gridsize = 11;
N = gridsize^2;
a0 = 0.5;
Rcell = 0.2*a0;
K = 10;
Con = 5;
flip = 1; % # cells to flip

% use hexagonal lattice
[dist, pos] = init_dist_hex(gridsize, gridsize);
dist_vec = a0*dist(1,:);
r = dist_vec(dist_vec>0); % exclude self influence
fN = sum(sinh(Rcell)*sum(exp(Rcell-r)./r)); % calculate signaling strength

% generate cell_type (0 case type 1, 1 case type 2)
cell_type = zeros(N,1); % all the same here

% initialize ON cells
%p = 0.68;
%iniON = round(p*N);
iniON = 58;
cells_ini = zeros(N,1);
cells_ini(randperm(N,iniON)) = 1;

fname_str = strrep(sprintf('N%d_n%d_Con_%d_K_%d_a0_%.2f_flip_%d', ...
    N, iniON, Con, K, a0, flip), '.', 'p');
folder = 'H:\My Documents\Multicellular automaton\figures\spin_flip_final\trajectories';
%% Initial simulation
% Variables
hin = figure(1);
t = 0;
I = [];
Non = [];
%mom = [];
cells_hist = {};

cells = cells_ini;
% store initial values
update_cell_figure(hin, pos, a0, cells, cell_type, t);
%k = waitforbuttonpress;
Non(end+1) = sum(cells);
I(end+1) = moranI(cells, a0*dist);
cells_hist{end+1} = cells;
[cells_out, changed, ~] = update_cells(cells, dist, Con, K, a0, Rcell);
while changed
    t = t+1;
    %k = waitforbuttonpress;
    update_cell_figure(hin, pos, a0, cells_out, cell_type, t);
    cells_hist{end+1} = cells_out;
    I(end+1) = moranI(cells_out, a0*dist);
    Non(end+1) = sum(cells_out);
    cells = cells_out;
    [cells_out, changed, ~] = update_cells(cells, dist, Con, K, a0, Rcell);
end
cells_eq = cells_out;
%% Simulate flipped configuration
% generate new config with flipped cells 
cells = cells_eq;

idx = randperm(N, flip);
cells(idx) = ~cells(idx);

hin = figure(2);
t2 = 0;
I2 = [];
Non2 = [];
%mom2 = [];
cells_hist2 = {};

% store initial values
update_cell_figure(hin, pos, a0, cells, cell_type, t2);
%k = waitforbuttonpress;
Non2(end+1) = sum(cells);
I2(end+1) = moranI(cells, a0*dist);
cells_hist2{end+1} = cells;
[cells_out, changed, ~] = update_cells(cells, dist, Con, K, a0, Rcell);
while changed
    t2 = t2+1;
    %k = waitforbuttonpress;
    update_cell_figure(hin, pos, a0, cells_out, cell_type, t2);
    cells_hist2{end+1} = cells_out;
    I2(end+1) = moranI(cells_out, a0*dist);
    Non2(end+1) = sum(cells_out);
    cells = cells_out;
    [cells_out, changed, ~] = update_cells(cells, dist, Con, K, a0, Rcell);
end

%% Plot Non vs time before & after flipping
h3 = figure(3);
hold on
plot(0:t, Non, 'b-', 'LineWidth', 1.5);
plot(t+(0:t2), Non2, 'r-', 'LineWidth', 1.5);

plot(0, Non(1), 'bo');
plot(t, Non(end), 'bx');

plot(t+0, Non2(1), 'ro');
plot(t+t2, Non2(end), 'rx');

xlabel('t');
ylabel('Non');
set(gca,'FontSize', 24)

qsave = 0;
if qsave
    v=1;
    fname = fullfile(folder,...
        strcat(fname_str,'-v', int2str(v), '_Non_vs_t') );
    while exist(strcat(fname, '.pdf')) == 2
        v=v+1;
        fname = fullfile(folder,...
            strcat(fname_str,'-v', int2str(v), '_Non_vs_t') );
    end
    save_figure_pdf(h3, 10, 8, fname);
    save_figure_eps(h3, 10, 8, fname);
end
%% Save trajectories if interesting

qsave = 0;
if qsave
    close all
    v=1;
    fname = fullfile(folder, 'data',...
        strcat(fname_str,'-v', int2str(v), '.mat'));
    while exist(fname, 'file') == 2
        v=v+1;
        fname = fullfile(folder, 'data',...
            strcat(fname_str,'-v', int2str(v), '.mat'));
    end
    save(fname);
end