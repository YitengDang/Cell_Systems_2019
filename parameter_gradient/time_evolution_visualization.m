% Time evolution of a system with visualization of the dynamics without
% noise showing the count of nearest neighbors that are ON
close all
clear all
warning off
set(0, 'defaulttextinterpreter', 'latex');

%%
% lattice parameters
gz = 16;
N = gz^2;
a0 = 0.5;
rcell = 0.2;
Rcell = rcell*a0;
% circuit parameters
Con = 12;
K = 16;
% initial conditions
p0 = 0.5;
iniON = round(p0*N);
%iniON = 1;

% Initialize parameters
%[pos,ex,ey] = init_cellpos_hex(gridsize,gridsize);
%dist = dist_mat(pos,gridsize,gridsize,ex,ey);
mcsteps = 0;
[pos, dist, ~, ~] = initial_cells_random_markov_periodic(gz, mcsteps, rcell);

dist_vec = a0*dist(1,:);
r = dist_vec(dist_vec>0); % exclude self influence
fN = sum(sinh(Rcell)*sum(exp(Rcell-r)./r)); % calculate signaling strength

% generate cell_type (0 case type 1, 1 case type 2)
cell_type = zeros(N,1); % all the same here

%% Plot phase diagram
M_int = 1; lambda12 = 1;
h = plot_phase_diagram(gz, a0, rcell, M_int, K, Con, lambda12);

%% K profile
Ax = 0.5;
Ay = 0.5;
Lx = 1; %default size
lambda_x = 1/2*Lx;
lambda_y = 1/2*sqrt(3)/2*Lx;

% sinusoidal wave
% K_func = @(x, y) K.*(1 + Ax.*sin(2.*pi.*x/lambda_x) + Ay.*sin(2.*pi.*y/lambda_y));
% K_all = K_func(pos(:, 1), pos(:, 2));

% square wave
K_all = K.*(1 + Ax.*square(pos(:, 1).*(2*pi/lambda_x)).* Ay.*square(pos(:, 2).*(2*pi/lambda_y)));

%{
figure;
x = 0:0.01:1;
y = 0:0.01:1;
plot(pos(:, 1), K_all, 'bo');
%}

%{
figure;
imagesc(reshape(K_cells_all, gz, gz))
colorbar;
set(gca, 'YDir', 'normal');
%}

%% Plot profile
h = figure;
% all sizes in units of pixels
Sx = 800; %ax.Position(3); %512;
Sy = (sqrt(3)/2*(gz-1)+2)/(3/2*(gz+1))*Sx;
a0_fig = Sx/(sqrt(N)+1); %Lx/(3/2*(gz+1));
Rcell_fig = rcell*a0_fig;
set(h, 'Position', [100 100 Sx Sy]);

% set image properties
set(gca, 'YTick', [], 'XTick', [], 'Color', [0.8 0.8 0.8]);
title(gca, 'Simulate dynamics', 'FontSize', 20);
Lx = 1;
d = 2*rcell*Lx/(sqrt(N)+1);
Ly = sqrt(3)/2*Lx;
xlim([-d Lx+d]);
ylim([-d Ly+d]);

% --plot cells--
hold on
% colours
c_all = K_all;
clr_k = zeros(N, 3); % black boundaries
%markers = {'o', 's'};

plot_handle = scatter(pos(:,1), pos(:,2), Rcell_fig^2, c_all, 'filled', 'o');
scatter(pos(:,1), pos(:,2), Rcell_fig^2, clr_k, 'o'); % plot cell boundaries

% Plot box outline
plot([0 Lx], [0 0], 'k--');
plot([0 Lx], [Ly Ly], 'k--');
plot([0 0], [0 Ly], 'k--');
plot([Lx Lx], [0 Ly], 'k--');

% colorbar
c = colorbar;
%c.Ticks = 0:0.2:1;
set(c, 'FontSize', 12);
c.Label.String = '$$K(i)$$';
ylabel(c, '$$K$$', 'Interpreter', 'latex', 'FontSize', 16)
map = 'parula';
colormap(map);
pause(0.5);

%%
% initialize ON cells
cells = zeros(N,1);
cells(randperm(N,iniON)) = 1;
%cells(1:2:N) = 1;

% variables
t = 0;
I = [];
Non = [];
h = [];
cells_hist = {};

% initialize figure
hin = figure;
disp_mol = 1;
showI = 0;
plot_handle = reset_cell_figure(hin, pos, rcell);
update_figure_periodic_scatter(plot_handle, cells, t, disp_mol, showI, a0, dist)

%update_cell_figure(hin, pos, a0, cells, cell_type, t);
% save vars and update cells
cells_hist{end+1} = cells;
Non(end+1) = sum(cells);
I(end+1) = moranI(cells, a0*dist);
[cells_out, changed, h(end+1)] = update_cells(cells, dist, Con, K_all, a0, Rcell);
while changed
    pause(0.5);
    t = t+1;
    %k = waitforbuttonpress;
    update_figure_periodic_scatter(plot_handle, cells_out, t, disp_mol, showI, a0, dist)
    cells_hist{end+1} = cells_out;
    I(end+1) = moranI(cells_out, a0*dist);
    Non(end+1) = sum(cells_out);
    cells = cells_out;
    [cells_out, changed, h(end+1)] = update_cells(cells, dist, Con, K_all, a0, Rcell);
end
%% Plot trajectory in p, I space
%{
figure();
hold on
plot(Non/N, I, 'r');
plot(Non(1)/N, I(1), 'ro');
plot(Non(end)/N, I(end), 'rx');
xlim([0 1]);
ylim([-0.05 1]);
%}
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