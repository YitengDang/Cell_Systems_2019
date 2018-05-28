% This file compares the simulation and the analytical formula for the
% entropy of the population. These are lines in the entropy map that
% compares the simulation with the analytical formula.

clear variable
close all
warning off

% Parameters that define the entropy map file to be read
date_str = '20160906'; % Date tag of the file
hour_str = '0418'; % hour tag of the file
a0 = 0.5;
grid = 15;

% Lines to watch (indices of vector K and Son from the entropy map file)
K_plot = [3, 8];
Son_plot = [5, 10];

% Parameters of the plot
title_str = 'a_0^c > a_0 = 0.5';
save = 1;
% ticks on the graph
tick_K = [1 5 10 15 20];
tick_Son = [1 5 10 15 20 25 30];
% The type of line to be plot, must have the same number of elements as the
% largest vector between K_plot and Son_plot
line_type = {'--', '-'}; % analytical formula is plot as lines
mark_type = {'*', 'o'}; % Simulation is plot as marks

% Name of the output file
fig_file = 'compare';

% Load data
name = sprintf('%s_%s_Son_K_map_a0%d_N%d_hexagonal', date_str, hour_str, round(10*a0),grid);
fname = fullfile(pwd,'data', 'entropy_map', ...
    strcat(name, '.mat'));
load(fname)

% Initialize the hexagonal grid
N = gridsize^2;
[dist, pos] = init_dist_hex(gridsize, gridsize);

% Calculate the signaling strength
dist_vec = dist(1,:);
r = a0*dist_vec(dist_vec>0); % exclude self influence
fN = sum(sinh(Rcell)*sum(exp(Rcell-r)./r)); % calculate signaling strength sphere

% Plot the entropy map
h = figure(1);
im = imagesc(K, Son, log(omega_sim));
c = colorbar;
colormap('summer')
hold on

% Insert the lines on the entropy map
for i = 1:numel(K_plot)
    plot([K(K_plot(i)) K(K_plot(i))],[1 tick_Son(end)],'--b', 'LineWidth', 2)
end
for i = 1:numel(Son_plot)
    plot([1 tick_K(end)],[Son(Son_plot(i)) Son(Son_plot(i))],'--r', 'LineWidth', 2)
end
hold off
% Adjust labels, ticks and fonts
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
    out_file = fullfile(pwd, 'figures', strcat(fig_file, name,'_map'));
    print(h, out_file,'-dpdf','-r0')
end

% Plot the comparison for fixed K
h1 = figure(2);
legend_str = {};
hold on
for i = 1:numel(K_plot)
    % calculate the analytical formula
    Son_vec = linspace(Son(1), Son(end), 100);
    % Plot analytical
    omega = zeros(size(Son_vec));
    for j = 1:numel(Son_vec)
        [omega(j),~] = entropy_eq_sphere(dist_vec, Son_vec(j), K(K_plot(i)), a0, Rcell);
    end
    plot(Son_vec, omega, strcat(line_type{i}, 'b'), 'LineWidth', 1.5)
    legend_str{end+1} = sprintf('Analytical - K = %.1f', K(K_plot(i)));
    % Plot simulation
    plot(Son, log(omega_sim(:, K == K(K_plot(i)))), ...
        strcat(mark_type{i}, 'b'), 'MarkerSize', 12)
    legend_str{end+1} = sprintf('Simulated - K = %.1f', K(K_plot(i)));
end
hold off
% Add legend and labels
legend(legend_str)
set(gca, 'FontSize', 20)
set(gca, 'xtick', tick_Son)
xlabel('C_{ON}', 'FontSize', 24)
ylabel('Entropy', 'FontSize', 24)

% Organize and save
if save > 0
    set(h1,'Units','Inches');
    set(h1, 'Position', [0 0 10 6 ])
    pos = get(h1,'Position');
    set(h1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    out_file = fullfile(pwd, 'figures', strcat(fig_file, name,'_fixedK'));
    print(h1, out_file,'-dpdf','-r0')
end

% Plot the comparison for fixed Son
h2 = figure(3);
legend_str = {};
hold on
for i = 1:numel(Son_plot)
    % calculate the analytical formula
    K_vec = linspace(K(1), K(end), 100);
    % Plot analytical
    omega = zeros(size(K_vec));
    for j = 1:numel(K_vec)
        [omega(j),~] = entropy_eq_sphere(dist_vec, Son(Son_plot(i)), K_vec(j), a0, Rcell);
    end
    plot(K_vec, omega, strcat(line_type{i}, 'r'), 'LineWidth', 1.5)
    legend_str{end+1} = sprintf('Analytical - C_{ON} = %.1f', Son(Son_plot(i)));
    % Plot simulation
    plot(K, log(omega_sim(Son == Son(Son_plot(i)),:)), ...
        strcat(mark_type{i}, 'r'), 'MarkerSize', 12)
    legend_str{end+1} = sprintf('Simulated - C_{ON} = %.1f', Son(Son_plot(i)));
end
hold off
% Add legend and labels
legend(legend_str, 'Location', 'northwest')
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
    out_file = fullfile(pwd, 'figures', strcat(fig_file, name, '_fixedCon'));
    print(h2, out_file,'-dpdf','-r0')
end