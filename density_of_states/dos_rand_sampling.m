%% Random lattice trials to find density of states or other parameters
%
% Average number of samples needed to sample all the states (Coupon
% collector's problem) ~ n log(n) ~ N*2^N
% For N=121, already N*2^N ~ 10^(38)
clear all;
close all;
%clc;

L = 11;
N = L^2;
a0 = 1.5;
[dist, ~] = init_dist_hex(L, L);

%% Run trials
ntrials = 10^6; %10^6;
Ivals = zeros(1,ntrials);
thetavals = zeros(1, ntrials);
pvals = zeros(1,ntrials);
for i=1:ntrials
    if mod(i, 10^4)==0
        disp(i);
    end
    cells = randi(2,N,1)-1; %random lattice
    %[Ivals(i), thetavals(i)] = moranI(cells, dist); % energy of lattice
    [Ivals(i), ~] = moranI(cells, a0*dist); % energy of lattice
    pvals(i) = sum(cells)/N;
    %disp(lattice);
end

% Study spacing between allowed values of I
%{
Ivals_sorted = sort(Ivals);
delI = Ivals_sorted(2:end) - Ivals_sorted(1:end-1);

figure();
histogram(delI);
xlim([0 0.00004]);

% 
figure();
plot(delI(100:end-100));
%}
%% Load data
load_folder = 'H:\My Documents\Multicellular automaton\figures\dos\random_sampling';
load_file = 'probIrandp_rand_N225_runs1000000_bins100_margin0p01';
load( fullfile(load_folder, load_file) );
%%
% Max. value of theta = 1 => theta normalized to 1 rather than fN
cells = zeros(N,1); 
[~, maxtheta] = moranI(cells, dist);

% Bin data
nbins = 100; %#bins for binning data AND for histogram
epsilon = 0.01; %bin margin around lowest and highest values
Iedges = linspace(min(Ivals)-epsilon, max(Ivals)+epsilon, nbins+1);
thetaedges = linspace(min(thetavals)-epsilon, max(thetavals)+epsilon, nbins+1);
[Ibins, Iedges] = histcounts(Ivals, Iedges);
[thetabins, thetaedges] = histcounts(thetavals, thetaedges);

%% Figures
set(0,'defaulttextinterpreter','tex');

% Histogram of P(I)
h01=figure();
histogram(Ivals, Iedges,'Normalization','probability');
set(gca, 'FontSize', 14);
xlabel('I');%, 'Interpreter', 'latex');
ylabel('$$P(I)$$');% , 'Interpreter', 'latex');

% Histogram of P(theta)
h02=figure();
histogram(thetavals, thetaedges,'Normalization','probability');
set(gca, 'FontSize', 14);
xlabel('$$\Theta$$');%, 'Interpreter', 'latex');
ylabel('$$P(\Theta)$$');% , 'Interpreter', 'latex');

%% 
% 2D histogram of P(p,I)
h1=figure();
histogram2(pvals, Ivals, 'Normalization', 'probability');
xlabel('p');
ylabel('I');
zlabel('$$P(p,I)$$'); %, 'Interpreter', 'latex');

%%
% 2D heat map
h2 = figure();
pedges = -1/2/N:1/N:1+1/2/N;
Iedges =  linspace(min(Ivals), max(Ivals), 20);
[pcounts, pedges, Iedges]  = histcounts2(pvals, Ivals, pedges, Iedges);
pcenters = (pedges(1:end-1)+pedges(2:end))/2;
Icenters = (Iedges(1:end-1)+Iedges(2:end))/2;
colormap('hot');
imagesc(pcenters, Icenters, pcounts'/ntrials);
c=colorbar;
c.Label.String = 'P(p,I)';
set(gca,'FontSize', 14);
set(gca,'YDir','normal')
xlabel('p');
ylabel('I');

%% 2D heat map log
h2b = figure();
pedges = -1/2/N:1/N:1+1/2/N;
Iedges = -0.1:0.01:0.15; %linspace(min(Ivals), max(Ivals), 20);
[pcounts, pedges, Iedges]  = histcounts2(pvals, Ivals, pedges, Iedges);
pcenters = (pedges(1:end-1)+pedges(2:end))/2;
Icenters = (Iedges(1:end-1)+Iedges(2:end))/2;
colormap('hot');
imagesc(pcenters, Icenters, log(pcounts'/ntrials));
c=colorbar;
c.Label.String = 'log(P(p,I))';
set(gca,'FontSize', 14);
xlabel('p');
ylabel('I');
set(gca, 'FontSize', 24, 'YDir', 'normal');
%  adjust axes 

%%
% Save figure and data
saveq = 1;
save_folder = 'H:\My Documents\Multicellular automaton\figures\dos\random_sampling';
if saveq
    fname = strrep(sprintf('probIrandp_rand_N%d_runs%d_bins%d',...
        L^2, ntrials, nbins), '.','p');
    % save histogram
    %out_fig = fullfile(pwd, 'figures', 'omega_pI', 'random_sampling', ...
    %    fname); 
    %savefig(h1, strcat(out_fig,'_hist.fig'));
    %save_figure_pdf(h1, 10, 6, strcat(out_fig,'_hist.pdf'));    
    % save heat map
    %out_fig2 = fullfile(pwd, 'figures', 'omega_pI', 'random_sampling', ...
    %    fname); 
    %savefig(h2, strcat(out_fig2,'_heat_map.fig'));
    %save_figure_pdf(h2, 10, 6, strcat(out_fig2, '_heat_map.pdf'));    
    % save heat map log scale
    out_fig2 = fullfile(save_folder, fname);
    savefig(h2b, strcat(out_fig2,'_heat_map_log.fig'));
    save_figure_pdf(h2b, 10, 6, strcat(out_fig2, '_heat_map_log.pdf')); 
    % Save data
    %out_file = fullfile(pwd, 'data', 'dos', 'random_sampling', ...
    %    strcat(fname,'.mat')); 
    out_file = fullfile(save_folder, strcat(fname,'.mat')); 
    save(out_file,'ntrials', 'L', 'nbins', 'Ivals', 'Ibins', 'pvals', 'pedges', 'Iedges');
end