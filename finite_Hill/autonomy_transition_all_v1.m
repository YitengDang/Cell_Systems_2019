% Compares trajectories of the finite and infinite Hill coefficient models
% with the same initial configuration
% Contains plots from (1) distance_metric_v2, (2) I2, (3) spread_cells
clear all
close all
set(0,'defaulttextinterpreter', 'latex');
%warning off
%%
% lattice parameters
gridsize = 11;
N = gridsize^2;
% circuit parameters
Con = 16;
K = 8;
hill = 2;
% initial conditions
initialID = 'binaryrand'; %'binary'

% Loop over a0
a0list_run = 5.4;
a0list = a0list_run; %5.25:0.005:5.3;

%%
% d variables
dmean = zeros(numel(a0list), 1);
dstd = dmean;
unif_frac = zeros(numel(a0list), 1);

% I2 variables
I2mean = zeros(numel(a0list), 1);
I2std = I2mean;
I2mean_ord = I2mean;
I2std_ord = I2mean;
I2ordcount = I2mean;

% spread variables
sigma_mean = zeros(numel(a0list), 2); % column 1: ON cells, column 2: OFF cells
sigma_std = sigma_mean;
p_mean = zeros(numel(a0list), 1);
p_std = p_mean;

% histograms
d_edges = 0:5:(floor(N/5)+1)*5;
d_hist = zeros(numel(a0list), numel(d_edges)-1);

I2_edges = -0.08:0.02:0.18;
I2_hist = zeros(numel(a0list), numel(I2_edges)-1);

% histogram
sigma_edges = 0:0.02:0.2;
sigma_hist_ON = zeros(numel(a0list), numel(sigma_edges)-1);
sigma_hist_OFF = sigma_hist_ON;

for idx = 1:numel(a0list)
 
% setup
a0 = a0list(idx);
Rcell = 0.2*a0;

% get fN
[dist, pos] = init_dist_hex(gridsize, gridsize);
dist = round(dist, 5);
dist_vec = a0*dist(1,:);
r = dist_vec(dist_vec>0); % exclude self influence
fN = sum(sinh(Rcell)*sum(exp(Rcell-r)./r)); % calculate signaling strength

% Get single cell / uniform lattice fixed points
fp = zeros(3, 1);
x0 = [0.03 0.2 0.65]; %estimates based on previous graph
%hfunc = @update_function_uniform;
hfunc2 = @(x) (((Con-1)*x + 1)*(1+fN))^hill/(K^hill+ (((Con-1)*x + 1)*(1+fN))^hill ) - x;
for i=1:3
    fp(i) = fzero(hfunc2, x0(i));
end

%% Load data, analyze and plot
% Path to search for the saved data. It searchs by the name, defined by the
% parameters chosen
path = fullfile(pwd, 'data', 'dynamics_transition 2017-10-25'); 
straux = '(\d+)';
straux2 = '(\w+)';

if ismember(a0, a0list_run)
    fpattern = strrep(sprintf('N%d_a0_%.2f_Con%.2f_K%.2f_hill%.2f_%s-v%s', ...
        N, a0, Con, K, hill, initialID, straux), '.', 'p');
else
    fpattern = strrep(sprintf('N%d_a0_%.3f_Con%.2f_K%.2f_hill%.2f_%s-v%s', ...
        N, a0, Con, K, hill, initialID, straux), '.', 'p');
end

% Get all file names in the directory
listing = dir(path);
num_files = numel(listing)-2; %first two entries are not useful
count = 0;
for i = 1:num_files
    filename = listing(i+2).name;
    % remove extension and do not include txt files
    [~,name,ext] = fileparts(filename);
    if strcmp(ext, '.mat')
        count = count + 1;
        names{count} = name;
    end
end

% Plot all final states
h1 = figure(1);
hold on

% Save Hamming distances between states
d_all = [];
I2fin = [];
unif_count = 0;
sigma_ON_all = [];
sigma_OFF_all = [];
p = []; %fraction ON cells
count = 0;
for i = 1:numel(names)
    [tokens, ~] = regexp(names{i},fpattern,'tokens','match');
    if numel(tokens) > 0
        disp(names{i}) % displays the file name
        count = count+1;
        % load the data
        load(fullfile(path,strcat(names{i},'.mat')));
        %cells_final = cells_hist{end};
        
        % d
        d = sum(abs(cells_ini - round(cells_final)));
        d_all(end+1) = d;
        
        % I2
        [~, theta] = moranI(cells_final, a0*dist);
        p_temp = sum(cells_final)/N;
        I2fin(end+1) = theta - (2*p_temp-1)^2;
        
        % count nonuniform
        if sum(round(cells_final))==0 || sum(round(cells_final))==N
            unif_count = unif_count + 1;
        end
        
        % get ON and OFF cells
        cells_ON = (cells_final > fp(2));
        p(end+1) = mean(cells_ON);
        if p(end) > 0 % only register if there are ON cells
            sigma_ON_all(end+1) = std(cells_final(cells_ON));
        end
        if p(end) < 1
            sigma_OFF_all(end+1) = std(cells_final(~cells_ON));
        end
        
        % Plot all cells
        figure(h1);
        plot(repmat(count, N, 1), cells_final, 'x');
    end
end
unif_frac(idx) = unif_count/count;

% store data
% d
d_hist(idx,:) = histcounts(d_all, d_edges);
dmean(idx) = mean(d_all);
dstd(idx) = std(d_all);

% I2
I2_hist(idx, :) = histcounts(I2fin, I2_edges); % Heat map of I2 values
I2mean(idx) = mean(I2fin); 
I2std(idx) = std(I2fin);

% I2: Ordered states: define as having sigma < eps
eps = 10^(-9);
I2mean_ord(idx) = mean(I2fin(abs(I2fin)>eps));
I2std_ord(idx) = std(I2fin(abs(I2fin)>eps));
I2ordcount(idx) = sum(abs(I2fin)>eps);

% spread
%[sigma_hist_ON(idx, :), sigma_edges] = histcounts(sigma_ON_all);
sigma_mean(idx, 1) = mean(sigma_ON_all);
sigma_mean(idx, 2) = mean(sigma_OFF_all);
sigma_std(idx, 1) = std(sigma_ON_all)/sqrt(numel(sigma_ON_all));
sigma_std(idx, 2) = std(sigma_OFF_all)/sqrt(numel(sigma_OFF_all));
p_mean(idx) = mean(p);
p_std(idx) = std(p)/sqrt(numel(p));

% Properties h1 (states of all cells)
h1;
hold on
plot([0 count+1], [fp(1) fp(1)], 'b--');
plot([0 count+1], [fp(2) fp(2)], 'r--');
plot([0 count+1], [fp(3) fp(3)], 'b--');
set(gca, 'FontSize', 24);
xlabel('Simulation number');
ylabel('Expression level');
xlim([0 count+1]);

qsave = 0;
if qsave
    fname_str = strrep(sprintf('all_cells_a0_%.3f_N%d_K%d_Con%d_hill_%.1f_%s_500runs',...
        a0, N, K, Con, hill, initialID), '.', 'p');
    fname = fullfile(pwd, 'figures', 'autonomy_transition', fname_str);
    save_figure_pdf(h1, 10, 8, fname);
    save_figure_eps(h1, 10, 8, fname);
end
close(h1);

%% histogram of Hamming distances
%{
h2 = figure(1+idx);
histogram(dall);
xlabel('$$d(X_{\infty}, X_{n=2})$$');
ylabel('$$P(d)$$');
title(sprintf('simulations: %d', numel(dall)));
set(gca,'FontSize', 24);
set(gcf, 'Units', 'Inches', 'Position', [0.5 0.5 10 8]);

qsave = 0;
if qsave
    fname_str = strrep(sprintf('Hamming_dist_ini_final_N%d_a0_%.2f_K%d_Con%d_hill_%.1f_%s',...
        N, a0, K, Con, hill, initialID), '.', 'p');
    fname = fullfile(pwd, 'figures', 'comparison_infinite_Hill', fname_str);
    save_figure_pdf(h3, 10, 8, fname);
    save_figure_eps(h3, 10, 8, fname);
end
%close(h1);
%}
end

%% Plot mean Hamming distance against a0
%{
h3=figure(3);
errorbar(a0list, dmean, dstd, 'LineWidth', 3);
xlabel('$$a_0$$');
ylabel('$$\langle d \rangle$$');
set(gca,'FontSize', 24);
set(gcf, 'Units', 'Inches', 'Position', [0.5 0.5 10 8]);
xlim([a0list(1)-0.05 a0list(end)+0.05]);
%ylim([0 70]);

qsave = 0;
if qsave
    fname_str = strrep(sprintf('Hamming_dist_avg_a0_%.1fto%.1f_N%d_K%d_Con%d_hill_%.1f_%s_500runs_v2',...
        a0list(1), a0list(end), N, K, Con, hill, initialID), '.', 'p');
    fname = fullfile(pwd, 'figures', 'autonomy_transition_Hamming_dist', fname_str);
    save_figure_pdf(h3, 10, 8, fname);
    save_figure_eps(h3, 10, 8, fname);
end
%% Plot fraction of uniform lattices against a0
h4=figure(4);
plot(a0list, unif_frac, 'o-', 'LineWidth', 3);
xlabel('$$a_0$$');
ylabel('Fraction uniform lattices');
set(gca,'FontSize', 24);
set(gcf, 'Units', 'Inches', 'Position', [0.5 0.5 10 8]);
xlim([a0list(1)-0.005 a0list(end)+0.005]);
%ylim([0 70]);

qsave = 0;
if qsave
    fname_str = strrep(sprintf('Frac_unif_lattices_a0_%.2fto%.2f_N%d_K%d_Con%d_hill_%.1f_%s_500runs_v2',...
        a0list(1), a0list(end), N, K, Con, hill, initialID), '.', 'p');
    fname = fullfile(pwd, 'figures', 'autonomy_transition_Hamming_dist', fname_str);
    save_figure_pdf(h4, 10, 8, fname);
    save_figure_eps(h4, 10, 8, fname);
end

%% Plot heat map of d_H
h5 = figure(5);
bincenters = (d_edges(1:end-1)+d_edges(2:end))/2;
nsim = numel(d_all); % number of simulations
imagesc(a0list, bincenters, d_hist'/nsim);
c = colorbar;
xlabel('$$a_0$$');
ylabel('$$d(X_{in}, X_{out})$$');
ylabel(c, 'Probability');
set(gca,'YDir', 'Normal', 'FontSize', 24);
xticks(5.0:0.1:5.9);
caxis([0 1]);

qsave = 0;
if qsave
    fname_str = strrep(sprintf('d_heat_map_a0_%.2fto%.2f_N%d_K%d_Con%d_hill_%.1f_%s_500runs',...
        a0list(1), a0list(end), N, K, Con, hill, initialID), '.', 'p');
    fname = fullfile(pwd, 'figures', 'autonomy_transition_Hamming_dist', fname_str);
    save_figure_pdf(h5, 12, 8, fname);
    save_figure_eps(h5, 12, 8, fname);
end


%% Plot mean I2 against a0
h6=figure(6);
errorbar(a0list, I2mean, I2std, 'LineWidth', 2);
xlabel('$$a_0$$');
ylabel('$$\langle I_2(t_f) \rangle$$');
set(gca,'FontSize', 24);
set(gcf, 'Units', 'Inches', 'Position', [0.5 0.5 10 8]);
%xlim([a0list(1)-0.05 a0list(end)+0.05]);
xlim([a0list(1)-0.005 a0list(end)+0.005]);
%ylim([0 70]);

qsave = 0;
if qsave
    fname_str = strrep(sprintf('I2_final_avg_a0_%.2fto%.2f_N%d_K%d_Con%d_hill_%.1f_%s_500runs',...
        a0list(1), a0list(end), N, K, Con, hill, initialID), '.', 'p');
    fname = fullfile(pwd, 'figures', 'autonomy_transition_I2_final', fname_str);
    save_figure_pdf(h6, 12, 8, fname);
    save_figure_eps(h6, 12, 8, fname);
end

%% Number of ordered final states
h7 = figure(7);
plot(a0list, I2ordcount, '-o');
xlabel('a0');
ylabel('Number of ordered final states');


%% Heat map
h5 = figure(8);
bincenters = (I2_edges(1:end-1)+I2_edges(2:end))/2;
nsim = numel(I2fin); % number of simulations
imagesc(a0list, bincenters, I2_hist'/nsim);
c = colorbar;
xlabel('$$a_0$$');
ylabel('$$I_2$$');
ylabel(c, 'Probability');
set(gca,'YDir', 'Normal', 'FontSize', 24);
caxis([0 1]);
xticks(5.0:0.1:5.9);

qsave = 0;
if qsave
    fname_str = strrep(sprintf('I2_final_heat_map_a0_%.2fto%.2f_N%d_K%d_Con%d_hill_%.1f_%s_500runs',...
        a0list(1), a0list(end), N, K, Con, hill, initialID), '.', 'p');
    fname = fullfile(pwd, 'figures', 'autonomy_transition_I2_final', fname_str);
    save_figure_pdf(h5, 12, 8, fname);
    save_figure_eps(h5, 12, 8, fname);
end

%% Plot fraction of ON cells
%
h9 = figure(9);
errorbar(a0list, p_mean, p_std, 'LineWidth', 3)
xlabel('$$a_0$$');
ylabel('Fraction ON cells');
set(gca,'FontSize', 22);
set(gcf, 'Units', 'Inches', 'Position', [0.5 0.5 10 8]);
xlim([a0list(1)-0.005 a0list(end)+0.005]);

qsave = 0;
if qsave
    fname_str = strrep(sprintf('Frac_ON_cells_avg_a0_%.1fto%.1f_N%d_K%d_Con%d_hill_%.1f_%s_500runs',...
        a0list(1), a0list(end), N, K, Con, hill, initialID), '.', 'p');
    fname = fullfile(pwd, 'figures', 'autonomy_transition', fname_str);
    save_figure_pdf(h9, 12, 8, fname);
    save_figure_eps(h9, 12, 8, fname);
end
%}
%% Load Mathematica data (mean field calculation)
dir='H:\My Documents\Multicellular automaton\Mathematica\fixed_points';
open(fullfile(dir, 'max_variance_cells_vs_a0_ON.xls'));
open(fullfile(dir, 'max_variance_cells_vs_a0_OFF.xls'));
%% Plot mean and std spread against a0
h10 = figure(10);
hold on
plot(a0list, sigma_mean(:,1), 'bo-', 'LineWidth', 3);
plot(a0list, sigma_mean(:,2), 'ro-', 'LineWidth', 3);
%errorbar(a0list, sigma_mean(:,1), sigma_std(:,1), 'b', 'LineWidth', 3); %ON cells
%errorbar(a0list, sigma_mean(:,2), sigma_std(:,2), 'r', 'LineWidth', 3); %OFF cells
xlabel('$$a_0$$');
ylabel('Spread $$\sigma(X)$$');
set(gca,'FontSize', 22);
set(gcf, 'Units', 'Inches', 'Position', [0.5 0.5 12 8]);
xlim([a0list(1)-0.05 a0list(end)+0.05]);
%ylim([0 70]);

% Plot loaded data
plot(a0list, sigmaON(1:end-2), 'bo--', 'LineWidth', 2.5);
plot(a0list, sigmaOFF(1:end-2), 'ro--', 'LineWidth', 2.5);

qsave = 1;
if qsave
    fname_str = strrep(sprintf('Spread_avg_a0_%.1fto%.1f_N%d_K%d_Con%d_hill_%.1f_%s_blueON_redOFF_500runs_v2',...
        a0list(1), a0list(end), N, K, Con, hill, initialID), '.', 'p');
    fname = fullfile(pwd, 'figures', 'autonomy_transition', fname_str);
    save_figure_pdf(h10, 12, 8, fname);
    save_figure_eps(h10, 12, 8, fname);
end
%}