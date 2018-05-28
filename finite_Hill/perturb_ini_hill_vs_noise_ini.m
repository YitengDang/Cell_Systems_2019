% Calculate the distance (L1 norm) between final states of two
% configurations differing only in initial configurations for a certain
% given noise strength
close all
clear variables
%warning off

% Parameters of the system
gridsize = 11;
N = gridsize^2;
a0 = 5.4;
Rcell = 0.2*a0;
K = 8;
Con = 16;
hill = 2;
prec = 8;
%noise = 0.01; 
noiselist = 10.^(-2:0.2:0);
nvar = numel(noiselist);
sampmethod = 'montecarlo';
p0 = 0.3;

% use hexagonal lattice
[dist, pos] = init_dist_hex(gridsize, gridsize);

% Calculate the signaling strength
dist_vec = dist(1,:);
r = a0*dist_vec(dist_vec>0); % exclude self influence
fN = sum(Rcell*sum(exp(Rcell-r)./r)); % calculate signaling strength
%kon = N*(K-fN-Con)/(Con-1)/fN;
%koff = N*(K-fN-1)/(Con-1)/fN;
p = (0:N)./N;

% filename for saving
fname_str = strrep(sprintf('vs_ini_noise_N%d_Con_%d_K_%d_a0_%.2f_p0_%.2f_noise10exp%.1f_%.1f', ...
    N, Con, K, a0, p0, log10(noiselist(1)), log10(noiselist(end)) ), '.', 'p');

%% calculate the map
%{
dH_prob = zeros(nvar, N+1);
dH_mean = zeros(nvar, 1);
dH_std = zeros(nvar, 1);
dI_mean = zeros(nvar, 1);
dI_std = zeros(nvar, 1);

for idx=1:numel(noiselist)
    noise = noiselist(idx);
    fprintf('Noise = %.4f \n', noise);
    [dH_count, dI, dH_mean(idx), dH_std(idx)] = ...
        spin_flip_dist_eq_parallel_fixed_pin_hill(dist, Con, K, a0, Rcell, p0, hill, prec, noise);
    dH_prob(idx, :) = dH_count/sum(dH_count);
    dI_mean(idx) = mean(mean(dI));
    dI_std(idx) = sum(std(dI, 0, 2).^2)/size(dI, 1);
end

%% save mat file
qsave = 0;
if qsave
    version = 1;
    fname = fullfile(pwd, 'figures', 'sensitivity_initial_cond_continuum', 'data',...
        strcat(fname_str,'-v', int2str(version), '.mat'));
    while exist(fname, 'file') == 2
        version=version+1;
        fname = fullfile(pwd, 'figures', 'sensitivity_initial_cond_continuum', 'data',...
            strcat(fname_str,'-v', int2str(version), '.mat'));
    end
    save(fname);
end
%}
%% Load saved data
path = fullfile(pwd, 'figures', 'sensitivity_initial_cond_continuum', 'data');
[fname, path, ~] = uigetfile(path);
load(fullfile(path,fname));
%% Hamming distance between final configurations
% between unflipped and flipped configurations
set(0, 'defaulttextinterpreter', 'latex');

% Plot average and std
h1 = figure(1);
p1=errorbar(noiselist, dH_mean/N, dH_std/N, 'o-', 'LineWidth', 2);
set(gca,'FontSize', 24)
set(get(p1,'Parent'), 'XScale', 'log');
xlabel('pertubation strength');
ylabel('$$d_H/N$$', 'Interpreter', 'latex');
set(gcf, 'Units', 'Inches', 'Position', [1 1 10 8]);
%}

% Save
save_fig = 0; % save figure? 0:no, 1: yes
if save_fig > 0
    fname = fullfile(pwd, 'figures', 'sensitivity_initial_cond_continuum', 'vs_noise',...
    	strcat(fname_str,'-v', int2str(version), '_dH_mean_std_errorbar'));
    save_figure_pdf(h1, 10, 8, fname);
    save_figure_eps(h1, 10, 8, fname);
end
%% Plot std only
h2 = figure(2);
p2 = plot(noiselist, dH_std/N, 'o-', 'LineWidth', 2);
set(gca,'FontSize', 24)
set(get(p2,'Parent'), 'XScale', 'log');
xlabel('pertubation strength');
ylabel('$$\sigma(d_H)$$', 'Interpreter', 'latex');
set(gcf, 'Units', 'Inches', 'Position', [1 1 10 8]);

% Save
save_fig = 0; % save figure? 0:no, 1: yes
if save_fig > 0
    %fname = fullfile(pwd, 'figures', 'single_spin_flip', 'vs_pin',...
    %	strcat(fname_str,'-v', int2str(version), '_Hamming_avg'));
    fname = fullfile(pwd, 'figures', 'sensitivity_initial_cond_continuum', 'vs_noise',...
    	strcat(fname_str,'-v', int2str(version), '_dH_std'));
    save_figure_pdf(h2, 10, 8, fname);
    save_figure_eps(h2, 10, 8, fname);
end
%% Plot separate distributions
h3 = figure(3);
hold on
for i=1:nvar
    plot(p, dH_prob(i, :));
end
set(gca,'FontSize', 24)
xlabel('$$d_H/N$$');
ylabel('$$P(d_H/N)$$');
legendCell = cellstr(num2str(log10(noiselist)', 'noise=%-.1f'));
legend(legendCell)
set(gcf, 'Units', 'Inches', 'Position', [1 1 10 8]);

%save 
save_fig = 0; % save figure? 0:no, 1: yes
if save_fig > 0
    fname = fullfile(pwd, 'figures', 'sensitivity_initial_cond_continuum', 'vs_noise',...
    	strcat(fname_str,'-v', int2str(version), '_dH_distributions'));
    save_figure_pdf(h3, 10, 8, fname);
    save_figure_eps(h3, 10, 8, fname);
end
%% Plot heat map
h4 = figure(4);
hold on
im_fig = imagesc(p, log10(noiselist), dH_prob);
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
ylabel('log(pertubation strength)', 'FontSize', 24)
xlim([0 0.6]);
ylim([log10(noiselist(1))-0.2 log10(noiselist(end))+0.2]);

%save 
save_fig = 0; % save figure? 0:no, 1: yes
if save_fig > 0
    fname = fullfile(pwd, 'figures', 'sensitivity_initial_cond_continuum', 'vs_noise',...
    	strcat(fname_str,'-v', int2str(version), '_dH_heat_map'));
    save_figure_pdf(h4, 10, 8, fname);
    save_figure_eps(h4, 10, 8, fname);
end
%% Return probability (dH=0)
h6 = figure(6);
p6 = plot(noiselist, dH_prob(:, 1), 'o-', 'LineWidth', 2);
xlabel('noise strength', 'FontSize', 24)
ylabel('$$P(d_H<1/N)$', 'FontSize', 24)
set(gca,'FontSize', 24)
set(get(p6,'Parent'), 'XScale', 'log');

%save 
save_fig = 0; % save figure? 0:no, 1: yes
if save_fig > 0
    fname = fullfile(pwd, 'figures', 'sensitivity_initial_cond_continuum', 'vs_noise',...
    	strcat(fname_str,'-v', int2str(version), '_prob_dH0'));
    save_figure_pdf(h6, 10, 8, fname);
    save_figure_eps(h6, 10, 8, fname);
end
%% Final I of trajectories
%
h5 = figure(5);
%plot(noiselist, dI_mean, 'r-o')
p5 = errorbar(noiselist, dI_mean, dI_std, 'r', 'LineWidth' ,2);
set(get(p5,'Parent'), 'XScale', 'log');

set(gca,'FontSize', 20)
%title(sprintf('$$N = %d, K = %.1f, S_{ON} = %.1f, a0 = %.1f, R = %.1f$$', ...
%    N, K, Con, a0, Rcell),'FontSize', 18, 'Interpreter', 'latex')
xlabel('pertubation strength', 'FontSize', 24)
ylabel('$$\langle dI \rangle$$', 'FontSize', 24)

%save
save_fig = 0; % save figure? 0: no, 1: yes
if save_fig > 0
    fname = fullfile(pwd, 'figures', 'sensitivity_initial_cond_continuum', 'vs_noise',...
    	strcat(fname_str,'-v', int2str(version), '_dI_mean_std'));
    save_figure_pdf(h5, 10, 8, fname);
    save_figure_eps(h5, 10, 8, fname);
end
%}