function plot_heat_map(X, Y, Z, fig_name, label_x, label_y, fig_title)
% Plot a heat map with all the required configurations
figure('Name', fig_name)
colormap('hot');
imagesc(X, Y, Z);
colorbar;
xlabel(label_x)
ylabel(label_y)
set(gca,'YDir','normal')
title(fig_title)