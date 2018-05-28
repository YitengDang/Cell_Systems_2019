% Analyse the plots of fraction of nonuniform lattices to get picture in
% K, Con space of the transition to autonomy
clear all 
close all
clc

% simulation parameters
N=121;
hill = 2;
Krange = 5:15;
Conrange = 10:25;
initialID = 'uniform';

% data
Klist_aut = [7 8 8 9 9 9 9 10 10 10 10 10 11 11 11 11 11 12 12 12 13];
Conlist_aut = [13 15 16 17 18 19 20 19 20 21 22 23 21 22 23 24 25 23 24 25 25];
a0_half = [5.25 4.5 5.5 4.25 4.75 5.25 6.25 4 4.5 4.75 5.25 5.75 3.75...
    4.25 4.5 4.75 5.25 3.75 4 4.25 3.75]; %value of a0 at which half of the lattices are non-uniform
a0_peak = [5 4.5 5.5 4.5 5 5.5 6 4 4.5 5 5.5 6 4 4.5 4.5 5 5.5 4 4 4.5 6]; % peak values of a0 for mean equilibration time
teq_vals = [145.05 150.65 366.08 128.68 106.18 255.48 94.48 112.13 148.01...
    101.60 81.98 86.54 119.71 498.96 79.15 106.51 126.44 147.72 83.26...
    74.92 83.17]; 
n = numel(Klist_aut);

% plot variables
Klist_aut_idx = zeros(1,n);
Conlist_aut_idx = zeros(1,n);
a0_map = zeros(numel(Conrange), numel(Krange));
teq_map = zeros(numel(Conrange), numel(Krange));

for i=1:n
    idxK = find(Klist_aut(i) == Krange);
    idxCon = find(Conlist_aut(i) == Conrange);
    Klist_aut_idx(i) = idxK;
    Conlist_aut_idx(i) = idxCon;
    a0_map(idxCon, idxK) = a0_half(i);
    teq_map(idxCon, idxK) = teq_vals(i);
end
%% Plot critical value of a0 
h=figure();
set(0, 'defaulttextinterpreter', 'latex');
imagesc(Krange, Conrange, a0_map);
c = colorbar;
colormap('summer')
xlabel('$$K$$', 'FontSize', 24)
ylabel('$$C_{ON}$$', 'FontSize', 24)
ylabel(c, '$$a_0^c$$', 'Interpreter', 'latex', 'FontSize', 24)
set(gca, 'ydir', 'normal', 'FontSize', 24)

qsave = 1;
if qsave
    fname_str = strrep(sprintf('Transition_a0c_N%d_K%.2fto%.2f_Con%.2fto%.2f_hill%.2f_%s',...
        N, Krange(1), Krange(end), Conrange(1), Conrange(end), hill, initialID), '.', 'p');
    out_file = fullfile(pwd, 'figures', 'finite_Hill_autonomy_KCon_map', fname_str);
    save_figure_pdf(h, 10, 8, out_file);
end
%% Plot maximum equilibration time
h=figure();
set(0, 'defaulttextinterpreter', 'latex');
imagesc(Krange, Conrange, teq_map);
c = colorbar;
colormap('hot')
set(gca,'FontSize', 24);
xlabel('$$K$$', 'FontSize', 24)
ylabel('$$C_{ON}$$', 'FontSize', 24)
ylabel(c, '$$max(t_{eq})$$', 'Interpreter', 'latex', 'FontSize', 24)
set(gca, 'ydir', 'normal', 'FontSize', 16)

qsave = 0;
if qsave
    fname_str = strrep(sprintf('Max_teq_N%d_K%.2fto%.2f_Con%.2fto%.2f_hill%.2f_%s',...
        N, Krange(1), Krange(end), Conrange(1), Conrange(end), hill, initialID), '.', 'p');
    out_file = fullfile(pwd, 'figures', 'finite_Hill_autonomy_KCon_map', fname_str);
    save_figure_pdf(h, 10, 8, out_file);
end