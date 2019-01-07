function h = plot_spatial_parameter(cData_all, parameter_str, title_str, rcell, pos)
% Plots profile of a parameter that varies in space as a heat map
%parameter_str = sprintf('$$K^{(%d %d)}$$', int_wave(1), int_wave(2));
N = size(cData_all, 1);
gz = sqrt(N);

h = figure;
% all sizes in units of pixels
Sx = 800; %ax.Position(3); %512;
Sy = (sqrt(3)/2*(gz-1)+2)/(3/2*(gz+1))*Sx;
a0_fig = Sx/(sqrt(N)+1); %Lx/(3/2*(gz+1));
Rcell_fig = rcell*a0_fig;
set(h, 'Position', [100 100 Sx Sy]);

% set image properties
set(gca, 'YTick', [], 'XTick', [], 'Color', [0.8 0.8 0.8]);
title(gca, title_str, 'FontSize', 20);
%title(gca, sprintf('$$A_x = %.1f, A_y = %.1f$$', Ax, Ay));
Lx = 1;
d = 2*rcell*Lx/(sqrt(N)+1);
Ly = sqrt(3)/2*Lx;
xlim([-d Lx+d]);
ylim([-d Ly+d]);

% --plot cells--
hold on
% colours
%c_all = K_all(:, int_wave(1), int_wave(2));
clr_k = zeros(N, 3); % black boundaries
%markers = {'o', 's'};

scatter(pos(:,1), pos(:,2), Rcell_fig^2, cData_all, 'filled', 'o');
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
%c.Label.String = '$$K(i)$$';
ylabel(c, parameter_str,...
    'Interpreter', 'latex', 'FontSize', 16)
map = 'parula';
colormap(map);

end