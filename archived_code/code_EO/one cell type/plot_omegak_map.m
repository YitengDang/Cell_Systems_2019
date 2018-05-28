% Plot the number of equilibrium states for different points in the
% phenotype map

clear variables
close all
warning off

% Parameters
date_str = '20160523';
hour_str = '1211';
title_str = 'a_0^c < a_0 = 1.5';
fig_file = 'phenotype_map_points_omegak_a0_15';
save = 1;
tick_K = [1 5 10 15 20];
tick_Son = [1 5 10 15 20 25 30];

% In which points would you like to plot omegak?
% point = (K, Son)

points = [3 24;6 21;10 21;15 20;17 14;20 14];
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
for i = 1:size(points,1)
    text(points(i,1),points(i,2),int2str(i), 'FontSize', 20)
end
K = 1:0.5:20;
plot([fN+1 fN+1], [1 Son(end)], 'b', 'LineWidth', 1.5) % ON region
plot(K, K./(1+fN), 'k', 'LineWidth', 1.5) % OFF region
plot(K, K-fN, 'g', 'LineWidth', 1.5) % autonomous 1
plot(K, (K-1)/fN, 'r', 'LineWidth', 1.5) % autonomous 2
hold off
%plot(K, K, '--k', 'LineWidth', 1.5) % Son = K
xlim([1 20])
ylim([1 30])
set(gca, 'ydir', 'normal', 'FontSize', 20)
set(gca, 'xtick', tick_K, 'ytick', tick_Son)
xlabel('K', 'FontSize', 24)
ylabel('C_{ON}', 'FontSize', 24)
ylabel(c, 'Entropy', 'FontSize', 24, 'rot', 90);
title(title_str, 'FontSize', 24)

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
for i = 1:size(points,1)
    % calculate the analytical formula
    [~,omegak] = entropy_eq_sphere(dist_vec, points(i,2), points(i,1), a0, Rcell);
    plot((0:N)/N , omegak, 'LineWidth', 1.5)
    legend_str{end+1} = int2str(i);
end
hold off
legend(legend_str, 'Location', 'eastoutside')
set(gca, 'FontSize', 20)
xlabel('p = n/N', 'FontSize', 24)
ylabel('p_{E, n}', 'FontSize', 24)

% Organize and save
if save > 0
    set(h1,'Units','Inches');
    set(h1, 'Position', [0 0 10 6 ])
    pos = get(h1,'Position');
    set(h1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    out_file = fullfile(pwd, 'figures', strcat(fig_file,'_pe'));
    print(h1, out_file,'-dpdf','-r0')
end