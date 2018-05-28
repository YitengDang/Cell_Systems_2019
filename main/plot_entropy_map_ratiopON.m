clear variables
close all

% Plot entropy map with the lines determining the autonomy region

% Parameters
date_str = '20160523';
hour_str = '0508';
title_str = 'a_0^c > a_0 = 0.5';
fig_file = 'entropy_map_ONcells_a005';
save = 1;
tick_K = [1 5 10 15 20];
tick_Son = [1 5 10 15 20 25 30];

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
[out, map] = get_phenotype_map(a0, dist, Rcell, Son, K);

% Plot the result
h = figure(1);
im = imagesc(K, Son, log(2^(N)*omega_sim));
c = colorbar;
colormap('summer')
hold on
Son= linspace(1,30,1000);
K = linspace(1,20,1000);
% pplus = 1;
% pminus = 0;
pplus = 0.5 + 1/sqrt(N)/2;
pminus = 0.5 - 1/sqrt(N)/2;
plot([fN+1 fN+1], [1 Son(end)], 'k', 'LineWidth', 1.5) % ON region
plot(K, K./(1+fN), 'k', 'LineWidth', 1.5) % OFF region
c1 = (K-fN*(1-pminus))/(1+pminus*fN);
c2 = (K-1-fN*(1-pplus))/pplus/fN;
if (c1(end/2) < c2(end/2))
    plot(K, c1, 'r', 'LineWidth', 1.5) % autonomous 1
    plot(K, c2, 'r', 'LineWidth', 1.5) % autonomous 2
end
%plot(K, K, '--k', 'LineWidth', 1.5) % Son = K
hold off
ylabel(c, 'Entropy', 'FontSize', 24, 'rot', 90);
ylim(c,[0, N*log(2)])
ylim([1 Son(end)])
xlim([1 K(end)])
title(title_str, 'FontSize', 24);
set(gca, 'ydir', 'normal', 'FontSize', 20)
set(gca, 'xtick', tick_K, 'ytick', tick_Son)
xlabel('K', 'FontSize', 24)
ylabel('C_{ON}', 'FontSize', 24)

figure(2)
him = imagesc(K, Son, out);
colormap(map);
set(gca, 'ydir', 'normal', 'FontSize', 16)

% Organize and save
if save > 0
    set(h,'Units','Inches');
    set(h, 'Position', [0 0 10 6 ])
    pos = get(h,'Position');
    set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    fig_file = fullfile(pwd, 'figures', fig_file);
    print(h, fig_file,'-dpdf','-r0')
end

