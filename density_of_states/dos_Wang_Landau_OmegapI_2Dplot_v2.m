% Plots results of Wang-Laudau algorithm for Omega(p,I) in 2D plots
% v2: takes into account different bin sizes
clear all;
close all;
clc;
set(0, 'defaulttextinterpreter', 'latex');
%% Parameters
L = 11; % size of lattice
N = L^2;
a0 = 1.5;
n = 12:61;

c1=1; f0 = exp(c1); % initial modification factor 
c2=4; ffinal = exp(10^(-c2)); % final modification factor
p_flat=0.8; % flatness parameter
qsave = 0;

% --Patch together results from different simulations --
splits = [12 24 36 48 55 61+1]; %points at which to break up n into smaller collections
nparts = numel(splits)-1; % number of parts
n_all = cell(1,nparts);
idx_all = cell(1,nparts); % indices corresponding to the n above, starting from 1 (instead of n(1))
nbins_sim = zeros(1, nparts);

for i=1:numel(splits)-1
    n_all{i} = splits(i):(splits(i+1)-1);
    idx_all{i} = n_all{i} - n(1) + 1;
    nbins_sim(i) = numel( n_all{i} ); % number of sims with given bin count
end
%n_all = {12:23, 24:35, 36:47, 48:54, 55:61};
%nbins_sim = [12 12 12 7 7]; % number of sims with given bin count

%% Load files
g_all = cell(1,nparts);
nbins = [36 36 39 39 39]; % number of histogram bins according to filename (might be wrong)
for i=1:nparts
    g_temp = zeros(nbins(i), nbins_sim(i));
    for j=1:nbins_sim(i)
        this_n = n(1) + sum(nbins_sim(1:i-1))+ (j-1);
        %fprintf('i=%d, j=%d', i,j);
        disp(this_n)
        fileid = strrep(sprintf('WL_norm_N%d_n0_%d_a0_%.2f_f0exp%d_ffin_e10e-%d_pflat%.1f_%dbins_fixed%.2f',...
            N, this_n, a0, c1, c2, p_flat, nbins(i)), '.', 'p');
        fname = fullfile('H:\My Documents\Multicellular automaton', 'data', 'dos', 'probI_fixedp', strcat(fileid, '.mat') ); % filename
        load(fname, 'lgn2');
        g_temp(:, j) = exp(lgn2)/sum(exp(lgn2));
    end
    g_all{i} = g_temp;
end
%% Plot P(I|p)
%{
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
omegapI = cell(1, numel(g_all));
for i=1:numel(g_all)
    idx = n_all{i} - n(1) + 1; % indices of omegap to take
    omegapI{i} = g_all{i}.*repmat(omegap(idx), nbins(i), 1);
end
%% Plot figures
h1=figure();
hold on
%colormap('hot');
%edge_max = [0.18 0.18 0.2 0.2 0.2];
edges = linspace(-0.06, 0.2, nbins(end)+1); % TUNE
Icenters = (edges(1:end-1)+edges(2:end))/2;
%Icenters = cell(1,numel(g_all));
omegapI_all = Inf*ones(max(nbins), numel(n)); %only works when smaller bins have same edges
for i=1:numel(g_all)-1
    %edges = linspace(-0.06, edge_max(i), nbins(i));
    %Icenters{i} = (edges(1:end-1)+edges(2:end))/2;
    omegapI_all(1:nbins(i), idx_all{i}) = omegapI{i};
    %imagesc(n_all{i}/N, Icenters{i}, log(omegapI{i})); 
    % (1) plot both sides
    if mod(N,2)
        plotdata = log([omegapI_all  fliplr(omegapI_all(:,1:end-1))] );
        hplot = imagesc([n/N (N+1-n(1:end-1))/N], Icenters, plotdata);
        %set(hplot, 'AlphaData', plotdata > 0); % doesn't work, tune range
    else
        plotdata = log([omegapI_all  fliplr(omegapI_all(:,1:end))]);
        hplot = imagesc([n/N (N+1-n)/N], Icenters, log([omegapI_all  fliplr(omegapI_all(:,1:end))]));
        %set(hplot, 'AlphaData', plotdata > 0);
    end
    % (2) Plot only one side
    %imagesc(n/N, Icenters, log(omegapI_all) ); 
    % Plot other side
end
c=colorbar;
c.Label.String = 'log(\Omega(p,I))';
set(gca,'FontSize', 24);
set(gca,'YDir','normal');
xlim([n(1)/N n(end)/N]);
ylim([-0.06 0.2]);
%xticks(0.1:0.1:0.9);
%yticks(edges(1:3:end));
xlabel('p');
ylabel('I');
%% Save data
if qsave
    fname = strrep(sprintf('OmegapI_WL_norm_N%d_a0_%.2f_f0exp%d_ffin_e10e-%d_pflat%.1f_max_%d_bins',...
        N, a0, c1, c2, p_flat, max(nbins)), '.', 'p');
    out_file = fullfile(pwd, 'data', 'dos', 'OmegapI_probI_fixedp', ...
        strcat(fname,'.mat')); % filename figure
    save(out_file);
end
%% Save heat map
if qsave
    fname = strrep(sprintf('OmegapI_WL_norm_full_N%d_a0_%.2f_f0exp%d_ffin_e10e-%d_pflat%.1f_max_%d_bins',...
        N, a0, c1, c2, p_flat, max(nbins)), '.', 'p');
    out_fig = fullfile(pwd, 'figures', 'dos', 'Wang-Landau probI_fixedp', ...
        strcat(fname,'.pdf')); % filename figure
    save_figure_pdf(h1, 10, 6, out_fig);
end