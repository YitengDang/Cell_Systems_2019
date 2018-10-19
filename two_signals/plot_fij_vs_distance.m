% Calculates f(r_ij), f_N and nearest neighbour approximations
% Plot results against parameters (N, a0, lambda)
close all
clear all
set(0,'defaulttextinterpreter', 'latex');
%% Parameters
% Fixed parameters
gz = 25;
N = gz^2;
a0 = 1.5;
rcell = 0.2;
Rcell = rcell*a0;
lambda = [1 1.2];

% get pos, dist
mcsteps = 0;
[pos, dist] = initial_cells_random_markov_periodic(gz, mcsteps, rcell);

% calculate fN
idx = gz + round(gz/2); % pick cell not at corner of grid
dist_vec = a0*dist(idx,:);
r = dist_vec(dist_vec>0);
fN = zeros(2,1);
fN(1) = sum(sinh(Rcell)*exp((Rcell-r)./lambda(1)).*(lambda(1)./r));
fN(2) = sum(sinh(Rcell)*exp((Rcell-r)./lambda(2)).*(lambda(2)./r));

% calculate f_ij for all unique distances
r_uniq = unique(dist_vec(dist_vec>0)); % exclude self influence
fij1 = sinh(Rcell)*exp((Rcell-r_uniq)./lambda(1)).*(lambda(1)./r_uniq) ; % calculate signaling strength
fij2 = sinh(Rcell)*exp((Rcell-r_uniq)./lambda(2)).*(lambda(2)./r_uniq) ; % calculate signaling strength

%% Plot f(rij) vs rij
h = figure;
hold on
plot(r_uniq/a0, fij1, 'bx-', 'LineWidth', 1.5);
plot(r_uniq/a0, fij2, 'rx-', 'LineWidth', 1.5);
xlabel('$r_{ij}/a_0$');
ylabel('$f(r_{ij})$');
set(gca, 'FontSize', 20);
xlim([1 ceil(max(r_uniq)/a0)]); 
set(h, 'Units', 'Inches', 'Position', [1 1 10 8]);

qsave = 0;
folder = 'H:\My Documents\Multicellular automaton\figures\two_signals\interaction_strength';
fname_str = strrep(sprintf('gz_%d_a0_%.2f_fij_vs_rij', gz, a0), '.', 'p');
fname = fullfile(folder,fname_str);
save_figure(h, 10, 8, fname, '.pdf', qsave)
%% Calculate nearest-neighbour approximation
nn_dist = [1 sqrt(3) 2]; % distance of the nearest neighbours
nn_count = [6 4 8]; % number of neighbours at given distance
r_uniq2 = nn_dist.*a0;
fij1 = sinh(Rcell)*exp((Rcell-r_uniq2)./lambda(1)).*(lambda(1)./r_uniq2) ; % calculate signaling strength
fij2 = sinh(Rcell)*exp((Rcell-r_uniq2)./lambda(2)).*(lambda(2)./r_uniq2) ; % calculate signaling strength

figure;
hold on
plot(r_uniq2/a0, fij1, 'bx-');
plot(r_uniq2/a0, fij2, 'rx-');

%% Plot NN interaction strength as bar graph
% Interaction strengths with nearest neighbour (first bar) and next to nearest neighbours (second bar)
h = figure;
hold on

fN_nn_first = [fij1(1)*nn_count(1); fij2(1)*nn_count(1)];
fN_nn_second = [sum(fij1(1:3).*nn_count(1:3)); sum(fij2(1:3).*nn_count(1:3)) ];
b1 = bar([0.8 1.2], fN_nn_first./fN, 'FaceColor','flat');
b2 = bar([1.8 2.2], fN_nn_second./fN, 'FaceColor','flat');

b1.CData(1,:) = [0 0 0.8];
b1.CData(2,:) = [0.8 0 0];
b2.CData(1,:) = [0 0 0.8];
b2.CData(2,:) = [0.8 0 0];

% custom legend
b1a = bar(NaN, 'FaceColor', [0 0 0.8]);
b2a = bar(NaN, 'FaceColor', [0.8 0 0]);
legend([b1a b2a], 'Signal 1', 'Signal 2');

set(gca, 'XTick', [1 2], 'XTickLabels', {'NN', 'NN+NNN'});
%plot([0 3], [fN(1) fN(1)], 'b--');
%plot([0 3], [fN(2) fN(2)], 'r--');
ylabel('$$f_N^{nn}/f_N$$');
set(gca, 'FontSize', 20);
set(h, 'Units', 'Inches', 'Position', [1 1 10 8]);
ylim([0 1]);

qsave = 0;
folder = 'H:\My Documents\Multicellular automaton\figures\two_signals\interaction_strength';
fname_str = strrep(sprintf('gz_%d_a0_%.2f_fij_nn_bar_plots', gz, a0), '.', 'p');
fname = fullfile(folder,fname_str);
save_figure(h, 10, 8, fname, '.pdf', qsave)

%% Plot fN vs parameter
gz = 25; a0 = 0.5; lambda = [1 1.2];
a0_all = [0.01 0.1:0.1:2];
gz_all = [5:5:35];
lambda_all = 1:0.1:2;
var_all = gz_all;

fN_all = zeros(numel(var_all), 2); % full
f_NN_all = zeros(numel(var_all), 2); % nearest neighbour
f_NNN_all = zeros(numel(var_all), 2); % next to nearest neighbour
for i=1:numel(var_all)
    %a0 = a0_all(i);
    gz = gz_all(i);
    %lambda(2) = lambda_all(i);
    
    % calculate pos, dist (only for varying N)
    mcsteps = 0;
    [pos, dist] = initial_cells_random_markov_periodic(gz, mcsteps, rcell);
    
    % calculate fN
    idx = gz + round(gz/2); % pick cell not at corner of grid
    dist_vec = a0*dist(idx,:);
    r = dist_vec(dist_vec>0);
    fN_all(i,1) = sum(sinh(Rcell)*exp((Rcell-r)./lambda(1)).*(lambda(1)./r));
    fN_all(i,2) = sum(sinh(Rcell)*exp((Rcell-r)./lambda(2)).*(lambda(2)./r));
    
    % calculate f_nn
    nn_dist = [1 sqrt(3) 2]; % distance of the nearest neighbours
    nn_count = [6 4 8]; % number of neighbours at given distance
    r_uniq2 = nn_dist.*a0;
    fij1 = sinh(Rcell)*exp((Rcell-r_uniq2)./lambda(1)).*(lambda(1)./r_uniq2) ; % calculate signaling strength
    fij2 = sinh(Rcell)*exp((Rcell-r_uniq2)./lambda(2)).*(lambda(2)./r_uniq2) ; % calculate signaling strength
    f_NN_all(i,:) = [fij1(1)*nn_count(1); fij2(1)*nn_count(1)];
    f_NNN_all(i,:) = [sum(fij1(1:3).*nn_count(1:3)); sum(fij2(1:3).*nn_count(1:3)) ];
end

h = figure;
hold on
plot(var_all, f_NN_all(:,1)./fN_all(:,1), 'bx-', 'LineWidth', 2);
plot(var_all, f_NN_all(:,2)./fN_all(:,2), 'rx-', 'LineWidth', 2);
plot(var_all, f_NNN_all(:,1)./fN_all(:,1), 'bx--', 'LineWidth', 2);
plot(var_all, f_NNN_all(:,2)./fN_all(:,2), 'rx--', 'LineWidth', 2);
%xlabel('$a_0$');
xlabel('Grid size ($\sqrt{N}$)');
%xlabel('$\lambda^{(2)}/\lambda^{(1)}$');
ylabel('Ratio');
%ylabel('$f_N^{(i)}$');
set(gca, 'FontSize', 20);
%xlim([1 ceil(max(r_uniq)/a0)]); 
ylim([0 1]);
legend({'f_{nn}^{(1)}/f_N^{(1)}', 'f_{nn}^{(2)}/f_N^{(2)}',...
    'f_{nnn}^{(1)}/f_N^{(1)}', 'f_{nnn}^{(2)}/f_N^{(2)}'},...
    'Location', 'ne');
set(gca, 'FontSize', 24);
set(h, 'Units', 'Inches', 'Position', [1 1 10 8]);

qsave = 1;
folder = 'H:\My Documents\Multicellular automaton\figures\two_signals\interaction_strength';
%fname_str = strrep(sprintf('plot_fnn_fN_ratio_vs_a0_gz%d_lambda12_%.1f', gz, lambda(2)), '.', 'p');
fname_str = strrep(sprintf('plot_fnn_fN_ratio_vs_gz_a0_%.1f_lambda12_%.1f', a0, lambda(2)), '.', 'p');
%fname_str = strrep(sprintf('plot_fnn_fN_ratio_vs_lambda12_gz%d_a0_%.1f', gz, a0), '.', 'p');
fname = fullfile(folder,fname_str);
save_figure(h, 10, 8, fname, '.pdf', qsave)