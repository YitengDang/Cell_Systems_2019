%% Random lattice trials to find density of states or other parameters
%
% Average number of samples needed to sample all the states (Coupon
% collector's problem) ~ n log(n) ~ N*2^N
% For N=121, already N*2^N ~ 10^(38)
clear all;
close all;
rng('shuffle');
clc;

L = 10;
N = L^2;
ntrials = 10^6; %trials for each p value
a0 = 1.5;

nbins = 100; %#bins for binning data AND for histogram
epsilon = 0.01; %bin margin around lowest and highest values

% Initialization
[dist, ~] = init_dist_hex(L, L);   
Ibincounts = zeros(N+1, nbins); 

% Single p value
%{
p=0.5;
iniON=round(p*N);
Ivals = zeros(ntrials,1);
for j=1:ntrials
    cells = zeros(L^2,1);
    cells(randperm(N, iniON))=1;
    Ivals(j) = moranI(cells, dist); % energy of lattice
    %pvals(i, j) = sum(cells)/N;
    %disp(lattice);
end
%}

% Loop over all p values
%
p = [50]/N; %(0:N)/N;
Ivals = zeros(numel(p), ntrials);
%pvals = zeros(numel(p),ntrials);
for i=1:numel(p)
    fprintf('p = %.2f \n', p(i) );
    % Run trials
    n = round( p(i)*N );
    for j=1:ntrials
        cells = zeros(L^2,1);
        cells( randperm(N, n) ) = 1;
        Ivals(i, j) = moranI(cells, dist); 
        %pvals(i, j) = sum(cells)/N;
        %disp(lattice);
    end
end
%}
%%
% Study spacing between allowed values of I
Ivals_sorted = sort(Ivals, 2);
delI = zeros(numel(p), size(Ivals, 2)-1);
for i=1:numel(p)
    delI(i, :) = Ivals_sorted(i, 2:end) - Ivals_sorted(i, 1:end-1);
end
%%
figure();
histogram(delI(2,:) );
%xlim([0 0.00004]);

%%
figure();
semilogy(1:size(delI,2), delI(1,:), 'x');
hold on
semilogy([1 size(delI,2)], [1/N 1/N], 'r--');



%%
% bin data
edges = linspace(min(min(Ivals))-epsilon, max(max(Ivals))+epsilon, nbins+1);
%edges = linspace(-0.06, 0.14, 31);
for i=1:N+1
    [Ibincounts(i, :), edges] = histcounts(Ivals(i, :), edges);
end

%%
% Plot heat map
set(0,'defaulttextinterpreter','latex');
h=figure();
colormap('hot');
imagesc(p(3:end-3), (edges(1:end-1)+edges(2:end))/2, Ibincounts(3:end-3, :)'/ntrials); 
%imagesc(p, (edges(1:end-1)+edges(2:end))/2, Ibincounts'/ntrials)
c=colorbar;
set(gca, 'FontSize', 14);
set(gca, 'YDir', 'normal');
c.Label.String = 'P(I|p)';
xlabel('p');
ylabel('I');
xticks(0.1:0.1:0.9);
yticks('auto');
%% Figures
% Histogram of P(I) averaged over p
figure();
histogram(Ivals, edges,'Normalization','probability');
set(gca, 'FontSize', 14);
xlabel('I');%, 'Interpreter', 'latex');
ylabel('$$P(I)$$');% , 'Interpreter', 'latex');
%%
% Save figure and data
saveq = 0;
if saveq
    fname = strrep(sprintf('probIfixedp_rand_fixedp_N%d_runs%d_bins%d_margin%.2f',...
        L^2, ntrials, nbins, epsilon), '.','p');
    % save heat map
    out_fig = fullfile(pwd, 'figures', 'dos', 'random_sampling', ...
        strcat(fname,'_heat_map2b.pdf'));
    save_figure_pdf(h, 10, 6, out_fig);
    
    % Save data
    close all;
    clear dist;
    out_file = fullfile(pwd, 'data', 'dos', 'random_sampling', ...
        strcat(fname,'.mat')); 
    save(out_file);
end 