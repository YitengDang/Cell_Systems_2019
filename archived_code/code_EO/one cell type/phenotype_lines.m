% Plots the limiting cases regions in color with the limiting lines

clear variables
close all

% Parameters of the system
gridsize = 21;
N = gridsize^2;
a0 = 1;
Rcell = 0.2*a0;
save_fig = 0;
fig_file = sprintf('sphere_a0_%d_N%d_regions', ceil(100*a0), N);

% use hexagonal lattice
if gridsize == 50
    load('dist5050.mat')
else
    [pos,ex,ey] = init_cellpos_hex(gridsize,gridsize);
    dist = dist_mat(pos,gridsize,gridsize,ex,ey);
end
dist_vec = dist(1,:);
r = a0*dist_vec(dist_vec>0); % exclude self influence

% Get the values for which to make the map
tick_K = [1 5 10 15 20];
tick_Son = [1 5 10 15 20 25 30];

% Get the signalling length
fN = sum(sinh(Rcell)*sum(exp(Rcell-r)./r)); % calculate signaling strength sphere

Son_vec = linspace(1,30,1000);
K_vec = linspace(1,20,1000);
[K, Son] = meshgrid(K_vec, Son_vec);

% Make 4 limiting regions as boolean matrices
R1 = (K-1 < fN); % Everything ON
R2 = (Son > K & (K-1)./Son > fN); % autonomous cells for Son > K
R3 = (K./Son - 1 > fN); % Everything OFF
R4 = (Son <= K & K - Son < fN & (K-1)./Son > fN); % autonomous cells for Son < K

h=figure(1);
out = R1 + 2*R2 + 3*R3 + 4*R4; % regions do not overlap
him = imagesc(K_vec, Son_vec, out);
set(him, 'AlphaData', out > 0);
% R1 -> red, R2 -> yellow, R3 -> magenta, R4 -> green, none -> white
map = [1, 1, 1
    1, 0, 0
    1, 1, 0
    1, 0, 1
    1, 1, 0];
tmp = map(1:max(max(out))+1,:);
colormap(tmp);

% Plot line indicating Son = K
hold on
Son= linspace(1,30,1000);
K = linspace(1,20,1000);
plot([fN+1 fN+1], [1 Son(end)], 'k', 'LineWidth', 1.5) % ON region
plot(K, K./(1+fN), 'b', 'LineWidth', 1.5) % OFF region
plot(K, K-fN, 'g', 'LineWidth', 1.5) % autonomous 1
plot(K, (K-1)/fN, 'r', 'LineWidth', 1.5) % autonomous 2
plot(K, K, '--k', 'LineWidth', 1.5) % Son = K
hold off
set(gca,'ydir','normal', 'FontSize', 24)
xlabel('K', 'FontSize', 30)
ylabel('C_{ON}', 'FontSize', 30)
ylim([1 Son(end)])
xlim([1 K(end)])
set(gca, 'xtick', tick_K, 'ytick', tick_Son)
title(sprintf('fN = %.3f, a_0 = %.2f', fN, a0), 'FontSize', 30)
box on

% Organize and save
if save_fig > 0
    set(h,'Units','Inches');
    set(h, 'Position', [0 0 9 6 ])
    pos = get(h,'Position');
    set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    fig_file = fullfile(pwd, 'figures', fig_file);
    print(h, fig_file,'-dpdf','-r0')
end

figure
[~,omegak] = entropy_eq_sphere(dist_vec, 5, 10, a0, Rcell);
for k = 0:N
    omegak(k+1) = omegak(k+1)*nchoosek(N, k);
end
plot((0:N)/N , omegak, 'LineWidth', 1.5)

figure
fa0 = sinh(Rcell)*sum(exp(Rcell-a0)./a0);
Kf = @(C, p, m, delta) 0.5.*(C+1) + fa0.*((C-1).*m + 1) + ...
    2.*(fN-6*fa0).*((C-1).*p+1) + delta/2;
plot(Kf(Son, 0.5, 3, 0), Son)
hold on
plot(Kf(Son, 1, 6, 0), Son)
plot(Kf(Son, 0, 0, 0), Son)
hold off
ylim([1 Son(end)])
xlim([1 K(end)])

diffY = @(C, K, p, m) K - (1 + fa0.*((C-1).*m + 1) + (fN-6*fa0).*((C-1).*p+1));
figure
plot(Son, 0.5*diffY(Son, Kf(Son, 0.3, 3, 0), 0.3, 3)./Kf(Son, 0.3, 3, 0))
