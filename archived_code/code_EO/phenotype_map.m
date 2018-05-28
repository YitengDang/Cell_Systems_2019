% Plots the limiting cases regions in color

clear variables
close all

% Parameters of the system
gridsize = 25;
N = gridsize^2;
a0 = 0.5;
Rcell = 0.2*a0;
save = 1;
fig_file = 'sphere_a0_05_N625_regions';

% use hexagonal lattice
[pos,ex,ey] = init_cellpos_hex(gridsize,gridsize);
dist = dist_mat(pos,gridsize,gridsize,ex,ey);
dist_vec = dist(1,:);
r = a0*dist_vec(dist_vec>0); % exclude self influence

% Get the values for which to make the map
Son_vec = linspace(1,30,1000);
K_vec = linspace(1,20,1000);
[K, Son] = meshgrid(K_vec, Son_vec);
tick_K = [1 5 10 15 20];
tick_Son = [1 5 10 15 20 25 30];

% Get the signalling length
fN = sum(Rcell*sum(exp(Rcell-r)./r)); % calculate signaling strength sphere
% fN = sum(sum(exp(-r)));

% Make 4 limiting regions as boolean matrices
R1 = (K-1 < fN); % Everything ON
R2 = (Son > K & (K-1)./Son > fN); % autonomous cells for Son > K
R3 = (K./Son - 1 > fN); % Everything OFF
R4 = (Son <= K & K - Son < fN & (K-1)./Son > fN); % autonomous cells for Son < K

h = figure(1);
out = R1 + 2*R2 + 3*R3 + 4*R4; % regions do not overlap
him = imagesc(K_vec, Son_vec, out);
%set(him, 'AlphaData', out > 0);
% R1 -> red, R2 -> yellow, R3 -> magenta, R4 -> green, none -> white
map = [1, 1, 1
    1, 0, 0
    1, 1, 0
    1, 0, 1
    0, 1, 0];
tmp = map(1:max(max(out))+1,:);
colormap(tmp);
%colorbar;

% Plot line indicating Son = K
hold on
plot(K_vec, K_vec, 'k--')
hold off
set(gca,'ydir','normal', 'FontSize', 16)
xlabel('K', 'FontSize', 18)
ylabel('S_{ON}', 'FontSize', 18)
set(gca, 'xtick', tick_K, 'ytick', tick_Son)

% Organize and save
if save > 0
    set(h,'Units','Inches');
    set(h, 'Position', [0 0 9 6 ])
    pos = get(h,'Position');
    set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    fig_file = fullfile(pwd, 'figures', fig_file);
    print(h, fig_file,'-dpdf','-r0')
end
