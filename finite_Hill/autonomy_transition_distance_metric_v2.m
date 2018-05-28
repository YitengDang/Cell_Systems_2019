% Compares trajectories of the finite and infinite Hill coefficient models
% with the same initial configuration
clear all
close all
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

%% Loop over a0
a0list = 5.25:0.005:5.3;
dmean = zeros(numel(a0list), 1);
dstd = dmean;
nonunif_frac = zeros(numel(a0list), 1);

% histogram
edges = 0:5:(floor(N/5)+1)*5;
d_hist = zeros(numel(a0list), numel(edges)-1);

for idx = 1:numel(a0list)
a0 = a0list(idx);
%% Load data, analyze and plot
% Path to search for the saved data. It searchs by the name, defined by the
% parameters chosen
path = fullfile(pwd, 'data', 'dynamics_transition 2017-10-25'); 
straux = '(\d+)';
straux2 = '(\w+)';

if a0 == 5.25 || a0 == 5.30
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
%h1 = figure(1);
%hold on

% Save Hamming distances between states
d_all = [];
nonunif_count = 0;
count = 0;
for i = 1:numel(names)
    [tokens, ~] = regexp(names{i},fpattern,'tokens','match');
    if numel(tokens) > 0
        disp(names{i}) % displays the file name
        count = count+1;
        % load the data
        load(fullfile(path,strcat(names{i},'.mat')));
        %cells_final = cells_hist{end};
        
        d = sum(abs(cells_ini - round(cells_final)));
        d_all(end+1) = d;
        
        %figure(h1);
        %plot(repmat(count, N, 1), cells_final, 'x');
        
        % count final states which are nonuniform
        if sum(round(cells_final))==0 || sum(round(cells_final))==N
            nonunif_count = nonunif_count + 1;
        end
    end
end
nonunif_frac(idx) = nonunif_count/count;

% store histogram counts
d_hist(idx,:) = histcounts(d_all, edges);
dmean(idx) = mean(d_all);
dstd(idx) = std(d_all);
%% histogram of Hamming distances
%{
set(0, 'defaulttextinterpreter', 'latex');
h3 = figure(1+idx);
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
set(0, 'defaulttextinterpreter', 'latex');
h100=figure(100);
errorbar(a0list, dmean, dstd, 'LineWidth', 2);
xlabel('$$a_0$$');
ylabel('$$\langle d \rangle$$');
set(gca,'FontSize', 24);
set(gcf, 'Units', 'Inches', 'Position', [0.5 0.5 10 8]);
xlim([a0list(1)-0.05 a0list(end)+0.05]);
%ylim([0 70]);

qsave = 0;
if qsave
    fname_str = strrep(sprintf('Hamming_dist_avg_a0_%.1fto%.1f_N%d_K%d_Con%d_hill_%.1f_%s',...
        a0list(1), a0list(end), N, K, Con, hill, initialID), '.', 'p');
    fname = fullfile(pwd, 'figures', 'autonomy_transition_Hamming_dist', fname_str);
    save_figure_pdf(h100, 10, 8, fname);
    save_figure_eps(h100, 10, 8, fname);
end
%% Plot fraction of uniform lattices against a0
h101=figure(101);
plot(a0list, nonunif_frac, 'o-', 'LineWidth', 2);
xlabel('$$a_0$$');
ylabel('Fraction uniform lattices');
set(gca,'FontSize', 24);
set(gcf, 'Units', 'Inches', 'Position', [0.5 0.5 10 8]);
xlim([a0list(1)-0.005 a0list(end)+0.005]);
%ylim([0 70]);

qsave = 0;
if qsave
    fname_str = strrep(sprintf('Frac_unif_lattices_a0_%.2fto%.2f_N%d_K%d_Con%d_hill_%.1f_%s',...
        a0list(1), a0list(end), N, K, Con, hill, initialID), '.', 'p');
    fname = fullfile(pwd, 'figures', 'autonomy_transition_Hamming_dist', fname_str);
    save_figure_pdf(h101, 10, 8, fname);
    save_figure_eps(h101, 10, 8, fname);
end

%% Plot heat map of d_H
h103 = figure(103);
bincenters = (edges(1:end-1)+edges(2:end))/2;
nsim = numel(d_all); % number of simulations
imagesc(a0list, bincenters, d_hist'/nsim);
c = colorbar;
xlabel('$$a_0$$');
ylabel('$$P(d(X_{in}, X_{out})|a_0)$$');
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