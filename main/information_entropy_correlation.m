% Plot correlation between multicellular entropy and mutual information
% (pin-pout)
clear variables
close all

% Load entropy values
%[fname, path, ~] = uigetfile(fullfile(pwd, 'data', 'entropy_map'));
path = fullfile(pwd, 'data', 'entropy_map');
fname = '20171114_0653_Son_K_map_a015_N15_hexagonal.mat';
load(fullfile(path, fname))
S_all = log(omega_sim);
%% Load mutual information values 
N = 225;
a0 = 1.5;
Klist = 1:30;
Conlist = 1:30;

fname_str = strrep(sprintf('Information_pin-pout_N%d_a0_%.1f_K_%dto%d_Con_%dto%d', ...
    N, a0, Klist(1), Klist(end), Conlist(1), Conlist(end)), '.', 'p');
fname = fullfile(pwd, 'data', 'pin_pout', 'analytical', strcat(fname_str, '.mat'));
load(fname);
%% scatter plot of data
h1 = figure(1);
plot(I_all, S_all, 'bo');
set(0, 'defaulttextinterpreter', 'latex');
xlabel('$$I(p_{out}, p_{in})$$');
ylabel('$$S = \log{(\Omega_E)}$$');
set(gca, 'FontSize', 24)

qsave = 0;
if qsave > 0 
    fname_str = strrep(sprintf('compare_Ip_S_N%d_a0_%.1f_K_%dto%d_Con_%dto%d_scatter',...
        N, a0, Klist(1), Klist(end), Conlist(1), Conlist(end)), '.', 'p');
    fname = fullfile(pwd, 'figures', 'information', 'compare_entropy', fname_str);
    save_figure_pdf(h1, 10, 8, fname);
    save_figure_eps(h1, 10, 8, fname);
end
%% Histogram I
figure();
histogram(I_all);
%% Histogram S
figure();
histogram(S_all);
%% 2D Histogram S and I
figure();
nbins = [10 10];
histogram2(I_all, S_all, nbins);
%% Heat map (joint probability density)
h2 = figure(2);
Xedges = linspace(0, 160, 11+1);
Yedges = linspace(0, 8, 8+1);
[counts, Sedges, Iedges] = histcounts2(S_all, I_all, Xedges, Yedges);
Ic = (Iedges(1:end-1)+Iedges(2:end))/2;
Sc = (Sedges(1:end-1)+Sedges(2:end))/2;
imagesc(Ic, Sc, counts/sum(sum(counts)) );
set(gca, 'YDir', 'normal', 'FontSize', 24);
c = colorbar;
%colormap('hot');
xticks(0:8);
xlabel('$$I(p_{out}, p_{in})$$');
ylabel('$$S = \log{(\Omega_E)}$$');
%ylabel(c, 'count');
ylabel(c, 'probability');

qsave = 0;
if qsave == 1 
    fname_str = strrep(sprintf('compare_Ip_S_N%d_a0_%.1f_K_%dto%d_Con_%dto%d_prob',...
        N, a0, Klist(1), Klist(end), Conlist(1), Conlist(end)), '.', 'p');
    fname = fullfile(pwd, 'figures', 'information', 'compare_entropy', fname_str);
    save_figure_pdf(h2, 10, 8, fname);
    save_figure_eps(h2, 10, 8, fname);
end