% Plots the limiting cases regions in color with the limiting lines. This
% overlays a color code (filling of the regions) with the colored lines.

clear variables
close all
set(0, 'defaulttextinterpreter', 'latex');

% Parameters of the system
gridsize = 15;
N = gridsize^2;
a0 = 1.5;
Rcell = 0.2*a0;
save_fig = 0;

% filename to print
fig_file = sprintf('sphere_a0_%d_N%d_regions', ceil(100*a0), N);

% Parameters to test
num_points = 1000;
Con_vec = linspace(1, 30, num_points);
K_vec = linspace(1, 20, num_points);
[K, Con] = meshgrid(K_vec, Con_vec);

% use hexagonal lattice
[dist, pos] = init_dist_hex(gridsize, gridsize);
dist_vec = dist(1,:);
r = a0*dist_vec(dist_vec>0); % exclude self influence

% Get the signalling length
fN = sum(sinh(Rcell)*sum(exp(Rcell-r)./r)); % calculate signaling strength sphere

% Make 4 limiting regions as boolean matrices
R1 = (1+fN - K) > 0; % Everything ON
R2 = ((1+fN)*Con - K ) < 0; % Everything OFF
R3 = ((Con + fN - K) > 0 & (1+fN - K) < 0); % ON remains ON & not all ON
R4 = ((1+ fN*Con - K) < 0 & ((1+fN)*Con - K ) > 0) ; % OFF remains OFF & not all OFF
%R3 = (Con > K & (K-1)./Con > fN); % autonomous cells for Son > K
%R4 = (Con <= K & K - Con < fN & (K-1)./Con > fN); % autonomous cells for Son < K

h=figure(1);
hold on
out = R1 + 2*R2 + 3*R3 + 4*R4; % only regions 3 and 4 can overlap
if ~isempty(find(unique(out)==0, 1))
    idx = 5; % activation-deactivation
    out(out==0) = 5; % ON remains ON & OFF remains OFF
elseif ~isempty(find(unique(out)==7, 1)) 
    idx = 6; % autonomy
    out(out==7) = 5; % ON remains ON & OFF remains OFF
end

him = imagesc(K_vec, Con_vec, out);
%set(him, 'AlphaData', out > 0); % invisible if not from any region
% R1 -> black
% R2 -> white
% R3 -> green
% R4 -> red
% activation-deactivation -> magenta
% autonomy -> gray
map = [0, 0, 0
    1, 1, 1
    0, 1, 0
    1, 0, 0
    1, 1, 0
    0.5, 0.5, 0.5];
tmp = map([1:4 idx], :);
%tmp = map(1:numel(unique(out)),:);
colormap(tmp);
%colormap(parula(5));
%}

% Plot lines
Con2 = linspace(1,30,1000);
K2 = linspace(1,20,1000);
%{
hold on
plot([fN+1 fN+1], [1 Con2(end)], 'k--', 'LineWidth', 1.5) % all ON region
plot([fN+1 fN+1], [1 Con2(end)], 'g--', 'LineWidth', 1.5) % all ON region
plot(K2, K2./(1+fN), 'b--', 'LineWidth', 1.5) % all OFF region
plot(K2, K2-fN, 'Color', [247 145 52]/256, 'LineStyle', '--', 'LineWidth', 1.5) % autonomous 1
plot(K2, (K2-1)/fN, 'r--', 'LineWidth', 1.5) % autonomous 2
% plot(K, 2*K/(1+fN)-1, '--k', 'LineWidth', 1.5) % (Con+1)(1+fN) = 2K
hold off
%}

% adjust the graph
set(gca,'ydir','normal', 'FontSize', 24)
xlabel('K', 'FontSize', 30)
ylabel('$$C_{ON}$$', 'FontSize', 30)
ylim([1 Con2(end)])
xlim([1 K2(end)])
title(sprintf('$$f_N = %.3f, a_0 = %.2f$$', fN, a0), 'FontSize', 30)

% set ticks
tick_K = [1 5 10 15 20];
tick_Son = [1 5 10 15 20 25 30];
set(gca, 'xtick', tick_K, 'ytick', tick_Son)

% Organize and save
if save_fig > 0
    set(h,'Units','Inches');
    set(h, 'Position', [0 0 9 6 ])
    pos = get(h,'Position');
    set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    fig_file = fullfile(pwd, 'figures', fig_file);
    print(h, fig_file,'-dpdf','-r0')
end
