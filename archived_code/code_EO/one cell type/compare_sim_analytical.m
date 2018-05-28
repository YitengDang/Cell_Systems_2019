% This file compares the simulation and the analytical formula for the
% entropy of the population. These are lines in the entropy map

clear variable
close all
warning off

% Parameters
date_str = '20160523';
hour_str = '1211';
title_str = 'a_0^c > a_0 = 1.5';
fig_file = 'sphere_a0_greater_a0c';
save = 1;
tick_K = [1 5 10 15 20];
tick_Son = [1 5 10 15 20 25 30];

K_plot = [5, 18];
Son_plot = [5, 21];
line_type = {'--', '-'};
mark_type = {'*', 'o'};

% Load data
fname = fullfile(pwd,'data',strcat(date_str,'_',hour_str,'_Son_K_map_N11_hexagonal.mat'));
load(fname)

% Get phenotype regions
N = gridsize^2;
[pos,ex,ey] = init_cellpos_hex(gridsize,gridsize);
dist = dist_mat(pos,gridsize,gridsize,ex,ey);

dist_vec = dist(1,:);
r = a0*dist_vec(dist_vec>0); % exclude self influence

% Get the signalling length
fN = sum(sinh(Rcell)*sum(exp(Rcell-r)./r)); % calculate signaling strength sphere

% Plot the result
h = figure(1);
im = imagesc(K, Son, log(2^(N)*omega_sim));
c = colorbar;
colormap('summer')
hold on
for i = 1:numel(K_plot)
    plot([K_plot(i) K_plot(i)],[1 Son(end)],'--b', 'LineWidth', 2)
end
for i = 1:numel(Son_plot)
    plot([1 K(end)],[Son_plot(i) Son_plot(i)],'--r', 'LineWidth', 2)
end
hold off
set(gca, 'ydir', 'normal', 'FontSize', 20)
set(gca, 'xtick', tick_K, 'ytick', tick_Son)
xlabel('K', 'FontSize', 24)
ylabel('C_{ON}', 'FontSize', 24)
ylabel(c, 'Entropy', 'FontSize', 24, 'rot', 90);

% Organize and save
if save > 0
    set(h,'Units','Inches');
    set(h, 'Position', [0 0 10 6 ])
    pos = get(h,'Position');
    set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    out_file = fullfile(pwd, 'figures', strcat(fig_file,'_map'));
    print(h, out_file,'-dpdf','-r0')
end

% K fixed
h1 = figure(2);
legend_str = {};
hold on
for i = 1:numel(K_plot)
    % calculate the analytical formula
    Son_vec = linspace(Son(1), Son(end), 100);
    omega = zeros(size(Son_vec));
    for j = 1:numel(Son_vec)
        [omega(j),~] = entropy_eq_sphere(dist_vec, Son_vec(j), K_plot(i), a0, Rcell);
    end
    plot(Son_vec, omega, strcat(line_type{i}, 'b'), 'LineWidth', 1.5)
    legend_str{end+1} = sprintf('Analytical - K = %d', K_plot(i));
    plot(Son, log(2^(N)*omega_sim(:, K == K_plot(i))), ...
        strcat(mark_type{i}, 'b'), 'MarkerSize', 12)
    legend_str{end+1} = sprintf('Simulated - K = %d', K_plot(i));
end
hold off
legend(legend_str)
set(gca, 'FontSize', 20)
set(gca, 'xtick', tick_K)
xlabel('C_{ON}', 'FontSize', 24)
ylabel('Entropy', 'FontSize', 24)

% Organize and save
if save > 0
    set(h1,'Units','Inches');
    set(h1, 'Position', [0 0 10 6 ])
    pos = get(h1,'Position');
    set(h1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    out_file = fullfile(pwd, 'figures', strcat(fig_file,'_fixedK'));
    print(h1, out_file,'-dpdf','-r0')
end

% Son fixed
h2 = figure(3);
legend_str = {};
hold on
for i = 1:numel(Son_plot)
    % calculate the analytical formula
    K_vec = linspace(K(1), K(end), 100);
    omega = zeros(size(K_vec));
    for j = 1:numel(K_vec)
        [omega(j),~] = entropy_eq_sphere(dist_vec, Son_plot(i), K_vec(j), a0, Rcell);
    end
    plot(K_vec, omega, strcat(line_type{i}, 'r'), 'LineWidth', 1.5)
    legend_str{end+1} = sprintf('Analytical - C_{ON} = %d', Son_plot(i));
    plot(K, log(2^(N)*omega_sim(Son == Son_plot(i),:)), ...
        strcat(mark_type{i}, 'r'), 'MarkerSize', 12)
    legend_str{end+1} = sprintf('Simulated - C_{ON} = %d', Son_plot(i));
end
hold off
legend(legend_str)
set(gca, 'FontSize', 20)
set(gca, 'xtick', tick_Son)
xlabel('K', 'FontSize', 24)
ylabel('Entropy', 'FontSize', 24)

% Organize and save
if save > 0
    set(h2,'Units','Inches');
    set(h2, 'Position', [0 0 10 6 ])
    pos = get(h2,'Position');
    set(h2,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    out_file = fullfile(pwd, 'figures', strcat(fig_file,'_fixedCon'));
    print(h2, out_file,'-dpdf','-r0')
end