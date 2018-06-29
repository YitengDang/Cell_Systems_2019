% Calculate the distance (L1 norm) between final states of two
% configurations differing only in initial configurations where 'flip'
% number of cells are flipped. 
close all
clear variables
%warning off

% Parameters of the system
gridsize = 11;
N = gridsize^2;
a0 = 0.5;
Rcell = 0.2*a0;
K = 10;
Con = 5;
flip_range = 1:10; 

% use hexagonal lattice
[dist, pos] = init_dist_hex(gridsize, gridsize);

% Calculate the signaling strength
dist_vec = dist(1,:);
r = a0*dist_vec(dist_vec>0); % exclude self influence
fN = sum(Rcell*sum(exp(Rcell-r)./r)); % calculate signaling strength
%kon = N*(K-fN-Con)/(Con-1)/fN;
%koff = N*(K-fN-1)/(Con-1)/fN;
p = (0:N)./N;

%% Load previous file to determine p_ini to consider
%
path = 'H:\My Documents\Multicellular automaton\figures\sensitivity_ini_spin_flip\data';
fname_in = fullfile(path, 'vs_pin_N121_Con_5_K_10_a0_5-v1.mat');
load(fname_in, 'dH_std');

figure();
plot(p, dH_std);
% find two largest elements in std of dH
[~, Nmax] = max(dH_std);
%}

% filename of the saved data
fname_str = strrep(sprintf('vs_flip_N%d_Con_%d_K_%d_a0_%.2f_n%d_flip_%dto%d', ...
    N, Con, K, a0, Nmax, flip_range(1), flip_range(end)), '.', 'p');

%% calculate the map
dH_prob = zeros(numel(flip_range), N+1);
dH_mean = zeros(numel(flip_range), 1);
dH_std = zeros(numel(flip_range), 1);
dI_mean = zeros(numel(flip_range), 1);
dI_std = zeros(numel(flip_range), 1);

for idx=1:numel(flip_range)
    flip = flip_range(idx);
    fprintf('Flip=%d \n', flip);
    kin = (Nmax(1)-1);
    [dH_count, dI, dH_mean(idx), dH_std(idx)] = ...
        spin_flip_dist_eq_parallel_fixed_pin(dist, Con, K, a0, Rcell, kin, flip);
    dH_prob(idx, :) = dH_count/sum(dH_count);
    dI_mean(idx) = mean(mean(dI));
    dI_std(idx) = sum(std(dI, 0, 2).^2)/size(dI, 1);
end

%% save mat file
qsave = 1;
if qsave
    version = 1;
    fname = fullfile(pwd, 'figures', 'sensitivity_ini_spin_flip', 'vs_flip', 'data',...
        strcat(fname_str,'-v', int2str(version), '.mat'));
    while exist(fname, 'file') == 2
        version=version+1;
        fname = fullfile(pwd, 'figures', 'sensitivity_ini_spin_flip', 'vs_flip', 'data',...
            strcat(fname_str,'-v', int2str(version), '.mat'));
    end
    save(fname);
end
%}
%% Hamming distance between final configurations
% between unflipped and flipped configurations
set(0, 'defaulttextinterpreter', 'latex');

% Plot average and std
h1 = figure(1);
errorbar(flip_range, dH_mean/N, dH_std/N, 'o-', 'LineWidth', 2);
set(gca,'FontSize', 24)
xlabel('flipped cells');
ylabel('$$d_H/N$$', 'Interpreter', 'latex');
set(gcf, 'Units', 'Inches', 'Position', [1 1 10 8]);
%}

% Save
save_fig = 1; % save figure? 0:no, 1: yes
if save_fig > 0
    fname = fullfile(pwd, 'figures', 'sensitivity_ini_spin_flip', 'vs_flip',...
    	strcat(fname_str,'-v', int2str(version), '_dH_mean_std_errorbar'));
    save_figure_pdf(h1, 10, 8, fname);
    save_figure_eps(h1, 10, 8, fname);
end
%% Plot std only
h2 = figure(2);
plot(flip_range, dH_std/N, 'o-', 'LineWidth', 2);
set(gca,'FontSize', 24)
xlabel(' flipped cells');
ylabel('$$\sigma(d_H)$$', 'Interpreter', 'latex');
set(gcf, 'Units', 'Inches', 'Position', [1 1 10 8]);

% Save
save_fig = 1; % save figure? 0:no, 1: yes
if save_fig > 0
    %fname = fullfile(pwd, 'figures', 'single_spin_flip', 'vs_pin',...
    %	strcat(fname_str,'-v', int2str(version), '_Hamming_avg'));
    fname = fullfile(pwd, 'figures', 'sensitivity_ini_spin_flip', 'vs_flip',...
    	strcat(fname_str,'-v', int2str(version), '_dH_std'));
    save_figure_pdf(h2, 10, 8, fname);
    save_figure_eps(h2, 10, 8, fname);
end
%% Plot separate distributions
h3 = figure(3);
hold on
for i=1:numel(flip_range)
    plot(p, dH_prob(i, :));
end
set(gca,'FontSize', 24)
xlabel('$$d_H/N$$');
ylabel('$$P(d_H/N)$$');
legendCell = cellstr(num2str(flip_range', 'flip=%-d'));
legend(legendCell)
set(gcf, 'Units', 'Inches', 'Position', [1 1 10 8]);

%save 
save_fig = 1; % save figure? 0:no, 1: yes
if save_fig > 0
    fname = fullfile(pwd, 'figures', 'sensitivity_ini_spin_flip', 'vs_flip',...
    	strcat(fname_str,'-v', int2str(version), '_dH_distributions'));
    save_figure_pdf(h3, 10, 8, fname);
    save_figure_eps(h3, 10, 8, fname);
end
%% Plot heat map
h4 = figure(4);
hold on
im_fig = imagesc(p, flip_range, dH_prob);
set(gca,'FontSize', 24)
% set title and font
%title(sprintf('$$N = %d, K = %.1f, S_{ON} = %.1f, a0 = %.1f, R = %.1f$$', ...
%    N, K, Con, a0, Rcell),'FontSize', 18, 'Interpreter', 'latex')
set(gca,'Ydir','normal','FontSize', 24)
% set invisible parts where count is zero
set(im_fig, 'AlphaData', dH_prob > 0);
% set colorbar and labels
c = colorbar;
c.Label.String = 'Probability';
xlabel('$$d_H/N$$', 'FontSize', 24)
ylabel('flipped cells', 'FontSize', 24)
xlim([0 0.6]);
ylim([flip_range(1)-0.5 flip_range(end)+0.5]);

%save 
save_fig = 1; % save figure? 0:no, 1: yes
if save_fig > 0
    fname = fullfile(pwd, 'figures', 'sensitivity_ini_spin_flip', 'vs_flip',...
    	strcat(fname_str,'-v', int2str(version), '_dH_heat_map'));
    save_figure_pdf(h4, 10, 8, fname);
    save_figure_eps(h4, 10, 8, fname);
end
%% Return probability (dH=0)
h6 = figure(6);
plot(flip_range, dH_prob(:, 1), 'o-', 'LineWidth', 2);
xlabel('flipped cells', 'FontSize', 24)
ylabel('$$P(d_H=0)$', 'FontSize', 24)
set(gca,'FontSize', 24)

%save 
save_fig = 1; % save figure? 0:no, 1: yes
if save_fig > 0
    fname = fullfile(pwd, 'figures', 'sensitivity_ini_spin_flip', 'vs_flip',...
    	strcat(fname_str,'-v', int2str(version), '_prob_dH0'));
    save_figure_pdf(h6, 10, 8, fname);
    save_figure_eps(h6, 10, 8, fname);
end
%% Final I of trajectories
%
h5 = figure(5);
%plot(flip_range, dI_mean, 'r-o')
errorbar(flip_range, dI_mean, dI_std, 'r')
set(gca,'FontSize', 20)
%title(sprintf('$$N = %d, K = %.1f, S_{ON} = %.1f, a0 = %.1f, R = %.1f$$', ...
%    N, K, Con, a0, Rcell),'FontSize', 18, 'Interpreter', 'latex')
xlabel('flipped cells', 'FontSize', 24)
ylabel('$$\langle dI \rangle$$', 'FontSize', 24)

%save
save_fig = 1; % save figure? 0: no, 1: yes
if save_fig > 0
    fname = fullfile(pwd, 'figures', 'sensitivity_ini_spin_flip', 'vs_flip',...
    	strcat(fname_str,'-v', int2str(version), '_dI_mean_std'));
    save_figure_pdf(h5, 10, 8, fname);
    save_figure_eps(h5, 10, 8, fname);
end
%}