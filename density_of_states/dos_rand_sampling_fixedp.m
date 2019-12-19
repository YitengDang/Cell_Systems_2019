%% Random lattice trials to find density of states or other parameters
%
% Average number of samples needed to sample all the states (Coupon
% collector's problem) ~ n log(n) ~ N*2^N
% For N=121, already N*2^N ~ 10^(38)
clear all;
close all;
rng('shuffle');
set(0,'defaulttextinterpreter','tex'); %works for all future plots
clc;
%%
L = 11;
N = L^2;
ntrials = 10^6; %trials for each p value
a0 = 1.5;

nbins = 100; %#bins for binning data AND for histogram
epsilon = 0.01; %bin margin around lowest and highest values

% Initialization
[dist, ~] = init_dist_hex(L, L);   

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
p_all = (0:N)/N;
Ivals = zeros(numel(p_all), ntrials);
%pvals = zeros(numel(p),ntrials);
for i=1:numel(p_all)
    fprintf('p = %.2f \n', p_all(i) );
    % Run trials
    n = round( p_all(i)*N );
    for j=1:ntrials
        cells = zeros(L^2,1);
        cells( randperm(N, n) ) = 1;
        Ivals(i, j) = moranI(cells, dist); 
        %pvals(i, j) = sum(cells)/N;
        %disp(lattice);
    end
end
%}
%% Load data
%
load_folder = 'M:\tnw\bn\hy\Shared\Yiteng\Multicellularity\data\main\dos\probI_fixedp_random_sampling';
load_file = 'probIfixedp_rand_fixedp_N121_runs1000000_bins100_margin0p01';
load( fullfile(load_folder, load_file) );
%}
%%
% Study spacing between allowed values of I
Ivals_sorted = sort(Ivals, 2);
delI = zeros(numel(p_all), size(Ivals, 2)-1);
for i=1:numel(p_all)
    delI(i, :) = Ivals_sorted(i, 2:end) - Ivals_sorted(i, 1:end-1);
end

figure();
histogram(delI(2,:) );
%xlim([0 0.00004]);
%{
figure();
semilogy(1:size(delI,2), delI(1,:), 'x');
hold on
semilogy([1 size(delI,2)], [1/N 1/N], 'r--');
%}

%%
% bin data
Ibincounts = zeros(numel(p_all), nbins); 
edges = linspace(min(min(Ivals))-epsilon, max(max(Ivals))+epsilon, nbins+1);
%edges = linspace(-0.06, 0.14, 31);
for i=1:numel(p_all) %N+1
    [Ibincounts(i, :), edges] = histcounts(Ivals(i, :), edges);
end

%%
% Plot heat map
h=figure();
Ibincenters = (edges(1:end-1)+edges(2:end))/2;

% all data
%{
prob = Ibincounts(:, :)'/ntrials;
image_h = imagesc(p_all, Ibincenters, prob); 
%xlim([0 1]);
%}

% partial data
%
p_offset_idx = 11; % index offset
prob = Ibincounts(p_offset_idx+1:end-p_offset_idx, :)'/ntrials;
image_h = imagesc(p_all(p_offset_idx+1:end-p_offset_idx-1), Ibincenters, prob); 
xlim([0.1 0.9]);
%}
%{
prob = Ibincounts(1:p_offset_idx, :)'/ntrials;
image_h = imagesc(p_all(1:p_offset_idx), Ibincenters, prob); 
xlim([0 0.1]);
%}
%image_h = imagesc(p_all(1:end), (edges(1:end-1)+edges(2:end))/2, Ibincounts(:, :)'/ntrials); 
set(image_h, 'AlphaData', prob>0);
%imagesc(p, (edges(1:end-1)+edges(2:end))/2, Ibincounts'/ntrials)
colormap('viridis');
%colormap('cool');
c=colorbar;
caxis([0 0.07]);
set(gca, 'FontSize', 24);
set(gca, 'YDir', 'normal');
c.Label.String = 'P(I|p)';
xlabel('p');
ylabel('I');
xticks(0:0.1:1);
yticks('auto');
%% Figures
%{
% Histogram of P(I) averaged over p
figure();
histogram(Ivals, edges,'Normalization','probability');
set(gca, 'FontSize', 14);
xlabel('I');%, 'Interpreter', 'latex');
ylabel('$$P(I)$$');% , 'Interpreter', 'latex');
%}
%%
% Save figure and data
saveq = 1;
save_folder = 'M:\tnw\bn\hy\Shared\Yiteng\Multicellularity\data\main\dos\probI_fixedp_random_sampling'; %'H:\My Documents\Multicellular automaton\figures\dos\random_sampling_fixedp';

if saveq
    fname_str = strrep(sprintf('probIfixedp_rand_fixedp_N%d_runs%d_bins%d_margin%.2f_p_0p1_to_0p9_viridis_v3_size_10_6',...
        L^2, ntrials, nbins, epsilon), '.','p');
    % save heat map
    %out_fig = fullfile(pwd, 'figures', 'dos', 'random_sampling', ...
    %    strcat(fname,'_heat_map2b.pdf'));
    out_fig = fullfile(save_folder, fname_str);
    save_figure_pdf(h, 10, 6, out_fig);
    
    % Save data
    close all;
    %{
    %clear dist;
    %out_file = fullfile(pwd, 'data', 'dos', 'random_sampling', ...
    %    strcat(fname_str, '.mat')); 
    out_file = fullfile(save_folder, strcat(fname_str, '.mat'));
    save(out_file, 'delI', 'edges', 'epsilon', 'Ibincounts', 'Ivals',...
        'Ivals_sorted', 'N', 'nbins', 'ntrials', 'p_all');
    %}
end 

%% ---(2) Plot Omega(p, I)---
% ---P(I|p) from random simulations, P(p) analytical---

% (a) Calculate P(p) analytically
% Calculate binomial coefficients (possibly not exact)
% first check if file exists already
fileid = sprintf('nchoosek_n%d', N);
fname = fullfile(pwd, 'data', 'nchoosek', strcat(fileid, '.mat'));
if exist(fname)==2
    %disp('exists!');
    load(fname);
else
    [out, msg] = mkdir('data', 'nchoosek');
    %disp('Doesnt exist!');
    omegap = zeros(1,N+1);
    for j=1:N+1  
        disp(j);
        % if not, calculate
        omegap(j) = nchoosek(N,j-1);
    end
    save(fname, 'omegap');
end
p_all = (0:N)/N;

% Plot Omega(p)
figure();
plot(p_all, omegap, 'x-');
xlabel('p');
ylabel('\Omega(p)'); %, 'Interpreter', 'latex');
%semilogy(p, omegap);
%% (b) Calculate Omega(p,I)
Icondp = Ibincounts/ntrials; % Conditional probability P(I|p), rows=p, columns=I 

% multiply each column of omegap with a row of Icondp to get P(I|p) P(p)
% for that value of p
omegapI = zeros(N+1, nbins);
for i=1:N+1
    omegapI(i,:) = Icondp(i,:)*omegap(i);
end

% (b2) set values below 1 equal to 0
idx = find(omegapI < 1);
omegapI(idx) = 0;
%%
% (c) Plot as surface
% mesh parameters
%nbins = 100; %#bins for binning data AND for histogram
%epsilon = 0.01;
Iedges = linspace(min(Ivals(:))-epsilon, max(Ivals(:))+epsilon, nbins+1);
Ibincenters = (Iedges(1:end-1)+Iedges(2:end))/2;
[X,Y] = meshgrid(p_all, Ibincenters);

%{
% plot surface
%set(0,'defaulttextinterpreter','latex');
h1=figure();
surf(X,Y,log(omegapI'));
%mesh(log(OmegapI));
xlabel('p');
ylabel('I');
zlabel('$$\log{\Omega(p,I)}$$'); %, 'Interpreter', 'latex');
%}

% (d) Plot as heat map
h2=figure();
colormap('viridis');
image_h = imagesc(p_all, Ibincenters, log(omegapI'));
c=colorbar;
c.Label.String = 'log(Omega(p,I))';
set(gca,'FontSize', 24);
set(image_h, 'AlphaData', omegapI'>0);
xlabel('p');
ylabel('I');
set(gca,'YDir','normal')
caxis([0 80]);

%% Save figures and data
% Save figure and data
saveq = 1;
save_folder = 'M:\tnw\bn\hy\Shared\Yiteng\Multicellularity\data\main\dos\OmegapI_probI_fixedp';
% 'H:\My Documents\Multicellular automaton\figures\dos\random_sampling_fixedp';
if saveq
    fname_str = strrep(sprintf('Omega_pI_rand_fixedp_N%d_runs%d_bins%d_margin%.2f_viridis_size_10_8',...
        L^2, ntrials, nbins, epsilon), '.','p');
    % save heat map
    %out_fig = fullfile(pwd, 'figures', 'dos', 'random_sampling', ...
    %    strcat(fname,'_heat_map2b.pdf'));
    out_fig = fullfile(save_folder, fname_str);
    save_figure_pdf(h2, 10, 8, out_fig);
    
    % Save data
    %{
    close all;
    %clear dist;
    %out_file = fullfile(pwd, 'data', 'dos', 'random_sampling', ...
    %    strcat(fname_str, '.mat')); 
    out_file = fullfile(save_folder, strcat(fname_str, '.mat'));
    save(out_file, 'edges', 'epsilon', 'Ibincounts', 'Ivals',...
        'Ivals_sorted', 'N', 'nbins', 'ntrials', 'p_all', 'omegap', 'omegapI');
    %}
end 