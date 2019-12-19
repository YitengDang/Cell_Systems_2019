% Plots results of Wang-Laudau algorithm for Omega(p,I) in 2D plots
clear all;
close all;
clc;
set(0, 'defaulttextinterpreter', 'tex');
%% Load files
L = 11; % size of lattice
N = L^2;
a0 = 0.5;
n = 11:61;
nbins = 30; % number of histogram bins
c1=1; f0 = exp(c1); % initial modification factor 
c2=4; ffinal = exp(10^(-c2)); % final modification factor
p_flat=0.8; % flatness parameter
qsave = 0;

g_all = zeros(nbins, numel(n));
%load_folder = fullfile(pwd, 'data', 'dos', 'probI_fixedp');
load_folder = 'M:\tnw\bn\hy\Shared\Yiteng\Multicellularity\data\main\dos\probI_fixedp\2017-06-01';
for i=1:numel(n)
    this_n = n(i);
    fileid = strrep(sprintf('WL_norm_N%d_n0_%d_a0_%.2f_f0exp%d_ffin_e10e-%d_pflat%.1f_%dbins_fixed%.2f',...
        N, this_n, a0, c1, c2, p_flat, nbins), '.', 'p');
    fname = fullfile(load_folder, strcat(fileid, '.mat') ); % filename
    load(fname, 'lgn2', 'edges');
    g_all(:, i) = exp(lgn2)/sum(exp(lgn2));
end
%% Plot P(I|p)
h2=figure;
Ibincenters = (edges(1:end-1)+edges(2:end))/2;
if mod(N,2)
    p_all = [n/N fliplr(N+1-n(1:end-1))/N];
    P_data = [g_all fliplr(g_all(:,1:end-1))];
else
    p_all = [n/N (N+1-n)/N];
    P_data = [g_all fliplr(g_all(:,1:end))];
end
colormap('viridis');
imagesc(p_all, Ibincenters, P_data);
c=colorbar;
c.Label.String = 'P(I|p)';
caxis([0 0.15]);
set(gca,'FontSize', 24);
set(gca,'YDir','normal');
xlabel('p');
ylabel('I');

% Save figure
qsave = 1;
save_folder = 'M:\tnw\bn\hy\Shared\Yiteng\Multicellularity\data\main\dos\probI_fixedp\2017-06-01';
if qsave
    fname_str = strrep(sprintf('probIfixedp_WL_norm_full_N%d_a0_%.2f_f0exp%d_ffin_e10e-%d_pflat%.1f_max_%d_bins_viridis_size_10_6',...
        N, a0, c1, c2, p_flat, max(nbins)), '.', 'p');
    fname = fullfile(save_folder, fname_str);
    save_figure_pdf(h2, 10, 6, fname);
end
%% Plot surface
set(0, 'defaulttextinterpreter', 'latex');
Ibincenters = (edges(1:end-1)+edges(2:end))/2;
[pm, Im] = meshgrid(n/N, Ibincenters);
h1=figure();
surf(pm, Im, g_all);
xlabel('p');
ylabel('I');

%% Plot heat map
h2=figure();
colormap('hot');
imagesc(n/N, Ibincenters, g_all);
c=colorbar;
c.Label.String = 'P(I|p)';
set(gca,'FontSize', 14);
set(gca,'YDir','normal');
xlabel('p');
ylabel('I');
%}
%% Calculate Omega(p,I)
% Get binomial coefficients (possibly not exact)
fileid = sprintf('nchoosek_n%d', N);
fname = fullfile(pwd,'data','nchoosek', strcat(fileid, '.mat'));
if exist(fname)==2
    %disp('exists!');
    load(fname);
else
    %disp('Doesnt exist!');
    binom = zeros(1,N+1);
    for j=1:N+1  
        disp(j);
        % if not, calculate
        binom(j) = nchoosek(N,j-1);
    end
    save(fname, 'binom');
end
%omegap = binom(n);
%% Plot figures
Icenters = (edges(1:end-1)+edges(2:end))/2;
h1=figure;
omegapI = g_all.*repmat(omegap(n), nbins, 1);
colormap('viridis');
X=sort(unique([n N-n(1:end-2)]))/N;
%patch values differently depending on whether N is even or odd
imagesc(X, Icenters, log([omegapI fliplr(omegapI(:, 1:end-mod(N,2)-1 ) )])); 
c=colorbar;
caxis([0 80]);
c.Label.String = 'log(\Omega(p,I))';
set(gca,'FontSize', 24);
set(gca,'YDir','normal');
%xticks(0.1:0.1:0.9);
%yticks(edges(1:3:end));
xlabel('p');
ylabel('I');

%% Save heat map
% Save figure
qsave = 1;
%save_folder = fullfile(pwd, 'figures', 'dos', 'Wang-Landau probI_fixedp');
save_folder = 'M:\tnw\bn\hy\Shared\Yiteng\Multicellularity\data\main\dos\probI_fixedp\2017-06-01';
if qsave
    fname = strrep(sprintf('OmegapI_WL_norm_N%d_a0_%.2f_f0exp%d_ffin_e10e-%d_pflat%.1f_%dbins_fixed%.2f_sizt_10_6',...
        N, a0, c1, c2, p_flat, nbins), '.', 'p');
    out_fig = fullfile(save_folder, strcat(fname,'.pdf')); % filename figure
    save_figure_pdf(h1, 10, 6, out_fig);
end