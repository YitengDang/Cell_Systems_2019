clear variables
close all
% Plot the entropy map for a saved data determined by the date and hour

% Parameters
save = 0;
tick_K = [1 5 10 15 20];
tick_Son = [1 5 10 15 20 25 30];

% Load data
[fname, path, ~] = uigetfile(fullfile(pwd,'data'));
load(fullfile(path,fname));

% Get phenotype regions
N = gridsize^2;
[dist, pos] = init_dist_hex(gridsize, gridsize);

dist_vec = dist(1,:);
r = a0*dist_vec(dist_vec>0); % exclude self influence

% Get the signalling length
fN = sum(sinh(Rcell)*sum(exp(Rcell-r)./r)); % calculate signaling strength sphere
%[out, map] = get_phenotype_map(a0, dist, Rcell, Son, K);

% Plot the result
h = figure(1);
im = imagesc(K, Son, omega_sim);
c = colorbar;
colormap('summer')
hold on
Son= linspace(1,30,1000);
K = linspace(1,20,1000);
plot([fN+1 fN+1], [1 Son(end)], 'b', 'LineWidth', 1.5) % ON region
plot(K, K./(1+fN), 'k', 'LineWidth', 1.5) % OFF region
%if (K(end/2)-fN) < (K(end/2)-1)/fN && (K(end/2)-1)/fN > Son(end/2)
    plot(K, K-fN, 'g', 'LineWidth', 1.5) % autonomous 1
    plot(K, (K-1)/fN, 'r', 'LineWidth', 1.5) % autonomous 2
%end
plot(K, K, '--k', 'LineWidth', 1.5) % Son = K
hold off
ylabel(c, 'Entropy = log(\Omega_E)', 'FontSize', 24, 'rot', 90);
ylim(c,[0, N*log(2)])
ylim([1 Son(end)])
xlim([1 K(end)])
%title(title_str, 'FontSize', 24);
set(gca, 'ydir', 'normal', 'FontSize', 20)
%set(gca, 'xtick', tick_K, 'ytick', tick_Son)
xlabel('K', 'FontSize', 24)
ylabel('C_{ON}', 'FontSize', 24)

% Organize and save
if save > 0
    set(h,'Units','Inches');
    set(h, 'Position', [0 0 10 6 ])
    pos = get(h,'Position');
    set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    fig_file = fullfile(pwd, 'figures', fig_file);
    print(h, fig_file,'-dpdf','-r0')
end

