%% Random lattice trials to find density of states or other parameters
%

clear all;
close all;
%clc;

L=11;
N = L^2;
[dist, ~] = init_dist_hex(L, L);

% Run trials
ntrials = 1000000;
Ivals = zeros(1,ntrials);
thetavals = zeros(1, ntrials);
pvals = zeros(1,ntrials);
for i=1:ntrials
    cells = randi(2,N,1)-1; %random lattice
    [Ivals(i), thetavals(i)] = moranI(cells, dist); % energy of lattice
    pvals(i) = sum(cells)/N;
    %disp(lattice);
end

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
set(0,'defaulttextinterpreter','latex');

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
[pcounts, pedges, Iedges]  = histcounts2(pvals, Ivals);
pbincenters = (pedges(1:end-1)+pedges(2:end))/2;
Ibincenters = (Iedges(1:end-1)+Iedges(2:end))/2;
colormap('hot');
imagesc(pbincenters, Ibincenters, pcounts'/ntrials);
c=colorbar;
c.Label.String = 'P(p,I)';
set(gca,'FontSize', 14);
xlabel('p');
ylabel('I');
%zlabel('$$\log{\Omega(p,I)}$$', 'Interpreter', 'latex');

%% 2D heat map log
h2b = figure();
[pcounts, pedges, Iedges]  = histcounts2(pvals, Ivals);
pbincenters = (pedges(1:end-1)+pedges(2:end))/2;
Ibincenters = (Iedges(1:end-1)+Iedges(2:end))/2;
colormap('hot');
imagesc(pbincenters, Ibincenters, log(pcounts'/ntrials));
c=colorbar;
c.Label.String = 'log(P(p,I))';
set(gca,'FontSize', 14);
xlabel('p');
ylabel('I');

%  adjust axes 

%%
% Save figure and data
saveq = 1;
if saveq
    fname = strrep(sprintf('OmegapI_rand_N%d_runs%d_bins%d_margin%.2f',...
        L^2, ntrials, nbins, epsilon), '.','p');
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
    out_fig2 = fullfile(pwd, 'figures', 'omega_pI', 'random_sampling', ...
        fname); 
    savefig(h2b, strcat(out_fig2,'_heat_map_log.fig'));
    save_figure_pdf(h2b, 10, 6, strcat(out_fig2, '_heat_map_log.pdf')); 
    % Save data
    out_file = fullfile(pwd, 'data', 'omega_pI', 'random_sampling', ...
        strcat(fname,'.mat')); 
    save(out_file,'ntrials', 'L', 'nbins', 'Ivals', 'Ibins', 'pvals');
end