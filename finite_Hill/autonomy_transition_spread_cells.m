% Compares trajectories of the finite and infinite Hill coefficient models
% with the same initial configuration
clear all
close all
%warning off
set(0, 'defaulttextinterpreter', 'latex');
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
a0list_run = 5.00:0.05:5.90;
a0list = a0list_run; %5.25:0.005:5.3;

sigma_mean = zeros(numel(a0list), 2); % column 1: ON cells, column 2: OFF cells
sigma_std = sigma_mean;
p_mean = zeros(numel(a0list), 1);
p_std = p_mean;

% histogram
edges = 0:0.02:0.2;
sigma_hist_ON = zeros(numel(a0list), numel(edges)-1);
sigma_hist_OFF = sigma_hist_ON;
%%
for idx = 1:numel(a0list)
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
        
        cells_ON = (cells_final > fp(2));
        p(end+1) = mean(cells_ON);
        if p(end) > 0 % only register if there are ON cells
            sigma_ON_all(end+1) = std(cells_final(cells_ON));
        end
        if p(end) < 1
            sigma_OFF_all(end+1) = std(cells_final(~cells_ON));
        end
        
        figure(h1);
        plot(repmat(count, N, 1), cells_final, 'x');
    end
end

% store histogram counts
%[sigma_hist_ON(idx, :), edges] = histcounts(sigma_ON_all);
sigma_mean(idx, 1) = mean(sigma_ON_all);
sigma_mean(idx, 2) = mean(sigma_OFF_all);
sigma_std(idx, 1) = std(sigma_ON_all);
sigma_std(idx, 2) = std(sigma_OFF_all);

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
%}
end

%% Plot fraction of ON cells
h100 = figure(100);
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
    save_figure_pdf(100, 12, 8, fname);
    save_figure_eps(100, 12, 8, fname);
end
%% Load Mathematica data (mean field calculation)
dir='H:\My Documents\Multicellular automaton\Mathematica\fixed_points';
open(fullfile(dir, 'max_variance_cells_vs_a0_ON.xls'));
open(fullfile(dir, 'max_variance_cells_vs_a0_OFF.xls'));
%% Plot mean and std spread against a0
h101 = figure(101);
hold on
errorbar(a0list, sigma_mean(:,1), sigma_std(:,1), 'b', 'LineWidth', 3); %ON cells
errorbar(a0list, sigma_mean(:,2), sigma_std(:,2), 'r', 'LineWidth', 3); %OFF cells
xlabel('$$a_0$$');
ylabel('Spread');
set(gca,'FontSize', 22);
set(gcf, 'Units', 'Inches', 'Position', [0.5 0.5 12 8]);
xlim([a0list(1)-0.05 a0list(end)+0.05]);
%ylim([0 70]);

% Plot loaded data
plot(a0list, sigmaON(1:end-2), 'bo--', 'LineWidth', 2.5);
plot(a0list, sigmaOFF(1:end-2), 'ro--', 'LineWidth', 2.5);

qsave = 0;
if qsave
    fname_str = strrep(sprintf('Spread_avg_a0_%.1fto%.1f_N%d_K%d_Con%d_hill_%.1f_%s_blueON_redOFF_500runs',...
        a0list(1), a0list(end), N, K, Con, hill, initialID), '.', 'p');
    fname = fullfile(pwd, 'figures', 'autonomy_transition', fname_str);
    save_figure_pdf(h101, 12, 8, fname);
    save_figure_eps(h101, 12, 8, fname);
end
%}
%% Plot heat maps of spread
% --- Fix normalization! ---------
%{
h102 = figure(102);
bincenters = (edges(1:end-1)+edges(2:end))/2;
nsim = size(sigma_hist_ON); % number of simulations
imagesc(a0list, bincenters, sigma_hist_ON'/nsim);
c = colorbar;
xlabel('$$a_0$$');
ylabel('$$P( \sigma |a_0)$$');
ylabel(c, 'Probability');
set(gca,'YDir', 'Normal', 'FontSize', 24);
caxis([0 1]);

h103 = figure(103);
bincenters = (edges(1:end-1)+edges(2:end))/2;
nsim = numel(sigma_all); % number of simulations
imagesc(a0list, bincenters, sigma_hist_OFF'/nsim);
c = colorbar;
xlabel('$$a_0$$');
ylabel('$$P( \sigma |a_0)$$');
ylabel(c, 'Probability');
set(gca,'YDir', 'Normal', 'FontSize', 24);
caxis([0 1]);

qsave = 0;
if qsave
    fname_str = strrep(sprintf('d_heat_map_a0_%.2fto%.2f_N%d_K%d_Con%d_hill_%.1f_%s',...
        a0list(1), a0list(end), N, K, Con, hill, initialID), '.', 'p');
    fname = fullfile(pwd, 'figures', 'autonomy_transition_Hamming_dist', fname_str);
    save_figure_pdf(h103, 14, 8, fname);
    save_figure_eps(h103, 14, 8, fname);
end
%}