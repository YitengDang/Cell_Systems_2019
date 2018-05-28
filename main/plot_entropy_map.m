clear variables
close all
% Plot the entropy map for a saved data determined by the date and hour

% Parameters
title_str = '$$a_0^c < a_0 = 1.5$$';
fig_file = 'sphere_a0_more_a0c';
tick_K = [1 5 10 15 20 25 30];
tick_Son = [1 5 10 15 20 25 30];

% Load data
% fname = fullfile(pwd,'data','entropy_map',strcat(date_str,'_',hour_str,'_Son_K_map_N15_hexagonal.mat'));
% load(fname)
[fname, path, ~] = uigetfile(fullfile(pwd,'data','entropy_map'));
load(fullfile(path, fname))

% Get phenotype regions
N = gridsize^2;
[dist, pos] = init_dist_hex(gridsize, gridsize);

dist_vec = dist(1,:);
r = a0*dist_vec(dist_vec>0); % exclude self influence

% Get the signalling length
fN = sum(sinh(Rcell)*sum(exp(Rcell-r)./r)); % calculate signaling strength sphere
[out, map] = get_phenotype_map(a0, dist, Rcell, Son, K);

% Plot the result
set(0, 'defaulttextinterpreter', 'latex');
h = figure(1);
im = imagesc(K, Son, log(omega_sim));
c = colorbar;
%colormap('summer')
hold on
Son = linspace(Son(1),Son(end),1000);
K = linspace(K(1),K(end),1000);
%plot([fN+1 fN+1], [1 Son(end)], 'b', 'LineWidth', 1.5) % ON region
%plot(K, K./(1+fN), 'k', 'LineWidth', 1.5) % OFF region
%if (K(end/2)-fN) < (K(end/2)-1)/fN && (K(end/2)-1)/fN > Son(end/2)
    %plot(K, K-fN, 'r', 'LineWidth', 1.5) % autonomous 1
    %plot(K, (K-1)/fN, 'LineWidth', 1.5, 'Color', [255 146 46]/256) % autonomous 2
%end
%plot(K, 2*K/(1+fN)-1, '--k', 'LineWidth', 1.5) % 0.5(Con_1)(1+fN) = K
hold off
ylabel(c, 'S = $$\log(\Omega_E)$$', 'FontSize', 24, 'rot', 90, 'Interpreter', 'latex');
ylim(c,[0, N*log(2)])
ylim([1 Son(end)])
xlim([1 K(end)])
title(title_str, 'FontSize', 24);
set(gca, 'ydir', 'normal', 'FontSize', 20)
%set(gca, 'xtick', tick_K, 'ytick', tick_Son)
xlabel('$$K$$', 'FontSize', 24)
ylabel('$$C_{ON}$$', 'FontSize', 24)

figure(2)
him = imagesc(K, Son, out);
colormap(map);
set(gca, 'ydir', 'normal', 'FontSize', 16)
%%
% Organize and save
save = 1;
if save > 0
    fname_str = strrep(sprintf('multicellular_entropy_map_K_Con_N%d_a0_%.1f', N, a0), '.', 'p');
    fname = fullfile(pwd, 'figures', 'information', fname_str); %filename
    save_figure_pdf(h, 10, 8, fname);
    save_figure_eps(h, 10, 8, fname);
end