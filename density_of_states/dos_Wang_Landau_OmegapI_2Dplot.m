% Plots results of Wang-Laudau algorithm for Omega(p,I) in 2D plots
clear all;
close all;
clc;
%% Load files
L = 11; % size of lattice
N = L^2;
a0 = 1.5;
n = 24:54;
nbins = 36; % number of histogram bins
c1=1; f0 = exp(c1); % initial modification factor 
c2=4; ffinal = exp(10^(-c2)); % final modification factor
p_flat=0.8; % flatness parameter
qsave = 0;

g_all = zeros(nbins, numel(n));
for i=1:numel(n)
    this_n = n(i);
    fileid = strrep(sprintf('WL_norm_N%d_n0_%d_a0_%.2f_f0exp%d_ffin_e10e-%d_pflat%.1f_%dbins_fixed%.2f',...
        N, this_n, a0, c1, c2, p_flat, nbins), '.', 'p');
    fname = fullfile(pwd, 'data', 'dos', 'probI_fixedp', strcat(fileid, '.mat') ); % filename
    load(fname, 'lgn2', 'edges');
    g_all(:, i) = exp(lgn2)/sum(exp(lgn2));
end
%% Plot P(I|p)
%
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

omegap = binom(n);
%% Plot figures
set(0, 'defaulttextinterpreter', 'latex');
Icenters = (edges(1:end-1)+edges(2:end))/2;
h1=figure();
omegapI = g_all.*repmat(omegap, nbins, 1);
colormap('hot');
X=sort(unique([n N-n(1:end-2)]))/N;
%patch values differently depending on whether N is even or odd
imagesc(X, Icenters, log([omegapI fliplr(omegapI(:, 1:end-mod(N,2)-1 ) )])); 
c=colorbar;
c.Label.String = 'log(\Omega(p,I))';
set(gca,'FontSize', 14);
set(gca,'YDir','normal');
xticks(0.1:0.1:0.9);
yticks(edges(1:3:end));
xlabel('p');
ylabel('I');
 
%% Save heat map
if qsave
    fname = strrep(sprintf('OmegapI_WL_norm_N%d_a0_%.2f_f0exp%d_ffin_e10e-%d_pflat%.1f_%dbins_fixed%.2f',...
        N, a0, c1, c2, p_flat, nbins), '.', 'p');
    out_fig = fullfile(pwd, 'figures', 'dos', 'Wang-Landau probI_fixedp', ...
        strcat(fname,'.pdf')); % filename figure
    save_figure_pdf(h1, 10, 6, out_fig);
end