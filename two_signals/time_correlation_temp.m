%% Computes the time correlation function for a trajectory
clear all
close all
clc
%% Load trajectory
folder = 'D:\data\videos_selected_data';
fname_str = 'chaos_limited_domain_to_death';
load(fullfile(folder, fname_str));

N = size(cells_hist{1}, 1);
gz = sqrt(N);
rcell = 0.2;
a0 = 1.5;
mcsteps = 0;

% Plot trajectory
%
[pos, dist] = initial_cells_random_markov_periodic(gz, mcsteps, rcell);
hin = figure;
plot_handle = reset_cell_figure(hin, pos, rcell);
t_ini = 0;
t_fin = numel(cells_hist)-1;
for t=t_ini:t_fin
    pause(0.01);
    cells = cells_hist{t+1};
    disp_mol = 12; showI = 1; 
    update_figure_periodic_scatter(plot_handle, cells, t, disp_mol, showI, a0, dist);
end
%}
%% Compute time correlation function
tmax = numel(cells_hist)-1;
l = 2;

%t1 = tmax;
corr_all = zeros(tmax+1, tmax+1, l, l); % correlation per interaction
corr_l = zeros(tmax+1, tmax+1); % full correlation

%i = 1; j =1; % correlation between i and j
for i=1:l^2
    [i1, i2] = ind2sub(2, i);
    for t1=0:tmax
        cells_1 = 2*cells_hist{t1+1}-1;
        for t=0:tmax
            cells = 2*cells_hist{t+1}-1;
            corr_all(t1+1, t+1, i1, i2) = sum(cells_1(:, i1).*cells(:, i2))/N;
                %- sum(cells_1(:, i1))/N*sum(cells(:, i2))/N;
            corr_l(t1+1, t+1) = sum(sum(cells_1.*cells))/N/l; %...
                %- sum(sum(cells_1))/N/l*sum(sum(cells))/N/l;
        end
    end
end
%% Plot for fixed t'
t1 = tmax;

h=figure;
colors = get(gca,'colororder');
hold on
%p1 = plot(0:tmax, corr_l, 'o-', 'LineWidth', 2);
p = zeros(1, l^2);
%corr_avg = zeros(tmax+1, 1);
for i=1:l^2
    [i1, i2] = ind2sub(2, i);
    if tmax<100
        plot(0:tmax, corr_all(t1+1, :, i1, i2), '--', ...
        'Color', colors(i, :), 'LineWidth', 1.5);
        p(i) = scatter(0:tmax, corr_all(t1+1, :, i1, i2), 'MarkerEdgeColor',...
            colors(i, :), 'MarkerFaceColor', ...
            colors(i, :), 'MarkerFaceAlpha', 0.2);
    else
        plot(0:tmax, corr_all(t1+1, :, i1, i2), '--', ...
        'Color', colors(i, :), 'LineWidth', 0.5);
        p(i) = scatter(0:tmax, corr_all(t1+1, :, i1, i2), 'MarkerEdgeColor',...
            colors(i, :), 'MarkerFaceColor', ...
            colors(i, :), 'MarkerFaceAlpha', 0.2);
    end
end
lw = (tmax<100) + 1;
corr_avg = squeeze(sum(sum(corr_all(t1+1,:,:,:), 4), 3))/l^2;
p2 = plot(0:tmax, corr_avg, 'k-', 'LineWidth', lw);
legend([p p2], {'11', '21', '12', '22', 'avg corr'});
%legend([p1 p], {'all', '11', '21', '12', '22'});
%% Heat map of correlations (avg)
t1_range = tmax-60:tmax; %0:tmax; %
t_range = tmax-60:tmax; %0:tmax; %

h=figure;
%corr_avg = squeeze(sum(sum(corr_all, 4), 3))/l^2;
corr_avg2 = (corr_all(:,:,1,1) + corr_all(:,:,2,2))./2;
imagesc(t1_range, t_range, corr_avg2(t1_range+1, t_range+1) );
xlabel('t');
ylabel('t''');
set(gca, 'YDir', 'normal', 'FontSize', 20);

% Add colorbar
c = colorbar;
c.Ticks = -1:0.2:1;
set(c, 'FontSize', 14);
ylabel(c, '$$\langle C^{(ij)}(t'', t) \rangle$$', 'Interpreter', 'latex',...
    'FontSize', 16, 'Rotation', 90);
%caxis([-1 1]);
%% Heat map of correlations (all)
% range to plot
t1_range = tmax-40:tmax;
t_range = tmax-40:tmax;

h = figure;
hold on
for i=1:l^2
    [i1, i2] = ind2sub(2, i);
    subplot(2,2,i);
    imagesc(t1_range, t_range, corr_all(1+t1_range, 1+t_range,i1, i2));
    xlabel('t');
    ylabel('t''');
    title(sprintf('C^{(%d %d)}(t'', t)', i1, i2));
    set(gca, 'Ydir', 'normal', 'FontSize', 14);
    caxis([-1 1]);
end
%{
% Add colorbar
c = colorbar;
c.Ticks = -1:0.2:1;
set(c, 'FontSize', 14);
ylabel(c, '$$\langle C^{(ij)}(t, t'') \rangle$$', 'Interpreter', 'latex',...
    'FontSize', 16, 'Rotation', 90);
caxis([-1 1]);
%}
%set(h, 'Units', 'Inches', 'Position', [0.1 0.1 11 6]);
set(h, 'units', 'normalized', 'outerposition', [0 0 1 1]);
