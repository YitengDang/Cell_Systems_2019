% This estimates based on the nearest neighbor approach:
% 1 - The minimum number of ON neighbors required to keep an ON cell ON
% 2 - The maximum number of ON cells possible to keep an OFF cell OFF

clear variables
close all

K = linspace(1, 20, 100);
m = 0:6;
% Parameters of the system
gridsize = 11;
N = gridsize^2;
a0 = 1.5;
Rcell = 0.2*a0;
save_fig = 1;

% use hexagonal lattice
[pos,ex,ey] = init_cellpos_hex(gridsize,gridsize);
dist = dist_mat(pos,gridsize,gridsize,ex,ey);
dist_vec = dist(1,:);
r = a0*dist_vec(dist_vec>0); % exclude self influence

% Get the signalling length
fN = sum(sinh(Rcell)*sum(exp(Rcell-r)./r)); % calculate signaling strength sphere
fa0 = sinh(Rcell)*sum(exp(Rcell-a0)./a0);

% Fractions to plot
p = [0.1 0.3 0.5 0.7 0.9];

[K_mesh, m_mesh] = meshgrid(K, m);

for i = 1:numel(p)
    aux = fa0*(m_mesh - 6*p(i));
    Con_ON = (K_mesh + aux - fN*(1-p(i)))./(aux + fN*p(i) + 1);
    Con_OFF = (K_mesh - 1 + aux - fN*(1-p(i)))./(aux + fN*p(i));
    h = figure(i);
    hold on
    plot(K, Con_ON, '--g', 'LineWidth', 1.5)
    plot(K, K-fN, 'g', 'LineWidth', 1.5)
    plot(K, Con_OFF, '--r', 'LineWidth', 1.5)
    plot(K, (K-1)/fN, 'r', 'LineWidth', 1.5)
    plot(K, K./(1+fN), 'b', 'LineWidth', 1.5)
    plot([fN+1 fN+1], [1 max(K)], 'k', 'LineWidth', 1.5)
    % plot(19, 14, 'ok') % plot a point
    hold off
    set(gca, 'FontSize', 20)
    xlabel('K', 'FontSize', 24)
    ylabel('C_{ON}', 'FontSize', 24)
    title(sprintf('p = %.2f', p(i)), 'FontSize', 24)
    xlim([1 max(K)])
    ylim([1 max(K)])
    box on
    if save_fig > 0
        filename = sprintf('N%d_a0%d_p%d', N, round(10*a0), round(100*p(i)));
        filename = fullfile(pwd, 'figures', 'nearest_neighbor', filename);
        save_figure_pdf(h, 8, 6, filename)
    end
end