clear all
close all
%% Load file
folder = fullfile('H:\My Documents\Multicellular automaton\app', 'data', 'time_evolution');
[file, path] = uigetfile(fullfile(folder, '\*.mat'), 'Load saved simulation');

if isnumeric(file) && isnumeric(path) % if user clicked "cancel"
    %app.MessagesTextArea.Value = sprintf('Simulation not loaded.');
    return
end
load(fullfile(path, file));
%% Plot single cell trajectories in X1-X2 space
t1 = 0; 
t2 = numel(cells_hist)-1;
cells_hist_red = cells_hist(t1+1:t2+1);

% convert data
tmax=numel(cells_hist_red);
cells_hist_mat_temp = cell2mat(cells_hist_red);
cells_hist_mat{1} = cells_hist_mat_temp(:, 1:2:end-1);
cells_hist_mat{2} = cells_hist_mat_temp(:, 2:2:end);

% cell to plot
N = size(cells_hist_red{1}, 1);
num_cells = 5; % choose num_cells to plot
cell_idx_all = randperm(N, num_cells);

for i=1:num_cells
    
cell_idx = cell_idx_all(1);

% plot cell trajectory (2D)
%{
figure;
hold on
plot(cells_hist_mat{1}(cell_idx,:), cells_hist_mat{2}(cell_idx,:), 'o--');
plot(cells_hist_mat{1}(cell_idx,1), cells_hist_mat{2}(cell_idx,1), 'ro', 'MarkerSize', 10);
plot(cells_hist_mat{1}(cell_idx,end), cells_hist_mat{2}(cell_idx,end), 'rx', 'MarkerSize', 10);
xlim([0 1]);
ylim([0 1]);
%}

%{
% plot cell trajectory (3D)
figure;
hold on
plot3(cells_hist_mat{1}(cell_idx,:), cells_hist_mat{2}(cell_idx,:), 0:tmax-1, 'o--');
plot3(cells_hist_mat{1}(cell_idx,1), cells_hist_mat{2}(cell_idx,1), 0, 'ro', 'MarkerSize', 10);
plot3(cells_hist_mat{1}(cell_idx,end), cells_hist_mat{2}(cell_idx,end), tmax-1, 'rx', 'MarkerSize', 10);
xlim([0 1]);
ylim([0 1]);
%}

% plot cell trajectory (3D) with arrows
figure;
hold on

plot3(cells_hist_mat{1}(cell_idx,:), cells_hist_mat{2}(cell_idx,:), 0:tmax-1, 'o--');
plot3(cells_hist_mat{1}(cell_idx,1), cells_hist_mat{2}(cell_idx,1), 0, 'ro', 'MarkerSize', 10);
plot3(cells_hist_mat{1}(cell_idx,end), cells_hist_mat{2}(cell_idx,end), tmax-1, 'rx', 'MarkerSize', 10);
xlim([-0.1 1.1]);
ylim([-0.1 1.1]);

% arrow endpoints at data points
%x = cells_hist_mat{1}(cell_idx,1:end-1);
%y = cells_hist_mat{2}(cell_idx,1:end-1);
%z = 0:tmax-2;
% arrow endpoints halfway lines
x = ( cells_hist_mat{1}(cell_idx,1:end-1) + cells_hist_mat{1}(cell_idx,2:end) )/2;
y = ( cells_hist_mat{2}(cell_idx,1:end-1) + cells_hist_mat{2}(cell_idx,2:end) )/2;
z = 0.5:tmax-1.5;
%plot3(x, y, z, 'o');

u = ( cells_hist_mat{1}(cell_idx,2:end) - cells_hist_mat{1}(cell_idx,1:end-1) );
v = ( cells_hist_mat{2}(cell_idx,2:end) - cells_hist_mat{2}(cell_idx,1:end-1) );
w = 0.5*ones(size(z));
quiver3(x, y, z, u, v, w, 0.01, 'b');
view(45, 45);

end