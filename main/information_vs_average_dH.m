% This script calculates mutual information and average dH (distance
% between initial and final states) for the same parameters and compares the two
close all
clear variables
warning off
set(0,'defaulttextinterpreter', 'latex');
%%
% Parameters of the system
gridsize = 15;
N = gridsize^2;
a0 = 1.5;
Rcell = 0.2*a0;
flip = 1; % number of cells to flip

% all parameters
K_list = [1:3];
Con_list = [2 2 2];

%--IMPORTANT---
final = 0;
% final = 0: sensitivity to initial conditions
% final = 1: stability of final conditions
%--------------

% use hexagonal lattice
[dist, pos] = init_dist_hex(gridsize, gridsize);

% Calculate the signaling strength
dist_vec = dist(1,:);
r = a0*dist_vec(dist_vec>0); % exclude self influence
fN = sum(Rcell*sum(exp(Rcell-r)./r)); % calculate signaling strength
%% (1) Load all files from directory
% output variables
dH_avg = [];
I_list = [];

% Load binomial distribution
fname_binom = fullfile('H:\My Documents\Multicellular automaton\data\nchoosek', strcat(num2str(N),'.mat'));
load(fname_binom);
binomprob = (nck/2^N);
clear nck;

% Get file names in the directory of dH files
subfolder = 'N225_a0_1p50_K_1to30_Con_1to30';
path = fullfile(pwd, 'figures', 'spin_flip_ini', subfolder, 'data');
listing = dir(path);
num_files = numel(listing)-2; %first two entries are not useful
count = 0;
for i = 1:num_files
    filename = listing(i+2).name;
    % remove extension and do not include txt files
    [~,name,ext] = fileparts(filename);
    if strcmp(ext, '.mat')
        count = count + 1;
        names{count} = name;
    end
end

Con_list = zeros(numel(names), 1);
K_list = Con_list;
for li=1:numel(names)
    disp(names{li});
    straux = '(\d+)';
    fpattern = strrep(sprintf('vs_pin_N%d_Con_%s_K_%s_a0_%.2f_flip_%d-v%s',...
        N, straux, straux, a0, flip, straux), '.', 'p');
    [tokens, ~] = regexp(names{li},fpattern,'tokens','match');
    Con = str2double(tokens{1}{1});
    K = str2double(tokens{1}{2});
    
    % Exclude the data points in the 'insulating phases'
    if 1+fN>K || (1+fN)*Con<K || (Con+fN>K && 1+fN*Con<K)
        disp('continue');
        continue;
    end
        
    Con_list(li) = Con;
    K_list(li) = K;

    % ---Load dH data------
    load(fullfile(path, names{li}), 'dHprob');

    % Compute average dH
    dH_all = (0:N);
    
    %dH_avg(end+1) = mean(dH_all*dHprob); % average over pin values uniformly
    dH_avg(end+1) = (N+1)*mean(dH_all'.*(dHprob*binomprob)); %average over binomial distribution
        
    % ----Load pin-pout data---------
    %N = 225;
    fname_str = strrep(sprintf('pin_pout_N%d_Con_%d_K_%d_a0_%d', N, Con, K, 10*a0), '.', 'p');
    folder = fullfile(pwd, 'figures', 'spin_flip_ini', 'data_pin-pout');
    load(fullfile(folder, fname_str), 'prob_pout', 'prob');
    
    % Plot pin-pout map
    %{
    figure(li);
    im_fig = imagesc((0:N)/N, (0:N)/N, prob);
    % set title and font
    title(sprintf('$$N = %d, K = %.1f, C_{ON} = %.1f, a_0 = %.1f, R = %.1f$$', ...
        N, K, Con, a0, Rcell),'FontSize', 18)
    set(gca,'Ydir','normal','FontSize', 20)
    % set invisible parts where count is zero
    set(im_fig, 'AlphaData', prob > 0);
    % set colorbar and labels
    c = colorbar;
    c.Label.String = 'Probability';
    xlabel('p_{in}', 'FontSize', 24)
    ylabel('p_{out}', 'FontSize', 24)
    %}
    % ---Calculate mutual information---
    % Compute S[ P(p_out) ]
    prob_pout = sum(prob, 2)/(N+1); % P(p_out)
    idx = prob_pout>0;
    S_pout = - sum(prob_pout(idx).*log2(prob_pout(idx)));

    % Compute S[ P(p_out | p_in) ]
    P_pin = 1/(N+1);
    S_pout_pin = zeros(N+1, 1);
    for i=1:N
        P_pout_pin = prob(:, i);
        idx2 = P_pout_pin > 0;
        S_pout_pin(i) = - sum(P_pout_pin(idx2).*log2(P_pout_pin(idx2)));
    end

    % Compute I(pin, pout)
    I_in_out = S_pout - sum(P_pin.*S_pout_pin);
    I_list(end+1) = I_in_out;
    %-----Close figures----------
    %close all
end
%% (2) Load files according to a previously defined list
%{
% output variables
dH_avg = zeros(numel(K_list), 1);
I_list = zeros(numel(K_list), 1);

for li=1:numel(K_list)
    %run=1;
    K = K_list(li);
    Con = Con_list(li);
    %% Load dH data
    %N = 121;
    fname_str = strrep(sprintf('vs_pin_N%d_Con_%d_K_%d_a0_%.2f_flip_%d-v%d',...      
        N, Con, K, a0, flip, 1), '.', 'p');
    folder = fullfile(pwd, 'figures', 'spin_flip_ini', 'N225_a0_1p50_K_1to30_Con_1to30', 'data');
    load(fullfile(folder, fname_str), 'dHprob');

    % Compute average dH
    dH_all = (0:N);
    dH_avg(li) = mean(dH_all*dHprob);

    %% Load pin-pout data
    %N = 225;
    fname_str = strrep(sprintf('pin_pout_N%d_Con_%d_K_%d_a0_%d', N, Con, K, 10*a0), '.', 'p');
    folder = fullfile(pwd, 'figures', 'spin_flip_ini', 'data_pin-pout');
    load(fullfile(folder, fname_str), 'prob_pout', 'prob');
    
    % ---Calculate mutual information---
    % Compute S[ P(p_out) ]
    prob_pout = sum(prob, 2)/(N+1); % P(p_out)
    idx = prob_pout>0;
    S_pout = - sum(prob_pout(idx).*log2(prob_pout(idx)));

    % Compute S[ P(p_out | p_in) ]
    P_pin = 1/(N+1);
    S_pout_pin = zeros(N+1, 1);
    for i=1:N
        P_pout_pin = prob(:, i);
        idx2 = P_pout_pin > 0;
        S_pout_pin(i) = - sum(P_pout_pin(idx2).*log2(P_pout_pin(idx2)));
    end

    % Compute I(pin, pout)
    I_in_out = S_pout - sum(P_pin.*S_pout_pin);
    I_list(li) = I_in_out;
end
close all;
%}
%% Plot average dH and I
h1=figure(1);
I_max = log2(N);
plot(dH_avg/N, I_list/I_max, 'bo', 'MarkerSize', 10);
xlabel('$$\langle \langle dH \rangle \rangle/N$$');
ylabel('$$I(p_{in}, p_{out})/I_{max}$$');
set(gca, 'FontSize', 24);

qsave = 1;
dataID = 'dataset1_binomprob'; % name data sets separately
if qsave
    fname_str = strrep(sprintf('%s_dH_avg_vs_information_N%d_a0_%.1f_Rcell_%.2f_flip_%d_normalized',...
        dataID, N, a0, Rcell, flip), '.', 'p');
    fname = fullfile(pwd, 'figures', 'spin_flip_ini', subfolder, fname_str); %filename
    save_figure_pdf(h1, 10, 8, fname);
    %save_figure_eps(h1, 10, 8, fname);
end
%% Plot histogram of <<dH>>
h2 = figure(2);
%edges = 0.45:0.01:0.55; %0.9:0.01:1.1; %0:0.1:1.5;
histogram(dH_avg/N, 'Normalization', 'probability');
xlabel('$$\langle \langle dH \rangle \rangle/N$$');
ylabel('Probability');
set(gca, 'FontSize', 24);
%xlim([0 1]);

qsave = 1;
if qsave
    fname_str = strrep(sprintf('%s_dH_histogram_N%d_a0_%.1f_Rcell_%.2f_flip_%d_normalized',...
        dataID, N, a0, Rcell, flip), '.', 'p');
    fname = fullfile(pwd, 'figures', 'spin_flip_ini', subfolder, fname_str); %filename
    %save_figure_eps(h2, 10, 8, fname);
    save_figure_pdf(h2, 10, 8, fname);
end
%% Plot location of K, Con on phase map
%{
h3 = figure(3);
hold on
Con_plot = 1:30;
K_plot = 1:30;
% Plot data
plot(K_list, Con_list, 'bx', 'MarkerSize', 10);
% Plot phase lines
plot([fN+1 fN+1], [Con_plot(1) Con_plot(end)+1], 'g', 'LineWidth', 1.5) % ON region
plot([K_plot K_plot(end)+1], [K_plot./(1+fN) (K_plot(end)+1)/(1+fN)], 'r', 'LineWidth', 1.5) % OFF region
if (K_plot(end/2)-fN) < (K_plot(end/2)-1)/fN && (K_plot(end/2)-1)/fN > Con_plot(end/2)
    plot([K_plot K_plot(end)+1], [K_plot K_plot(end)+1]-fN, 'k', 'LineWidth', 1.5) % autonomous 1
    plot([K_plot K_plot(end)+1], ([K_plot K_plot(end)+1]-1)/fN, 'k', 'LineWidth', 1.5) % autonomous 2
end
xlim([1 30]);
ylim([1 30]);
xlabel('$$K$$');
ylabel('$$C_{ON}$$');
set(gca, 'FontSize', 24);

qsave = 0;
if qsave
    fname_str = strrep(sprintf('%s_data_points_N%d_a0_%.1f_Rcell_%.2f_flip_%d',...
        dataID, N, a0, Rcell, flip), '.', 'p');
    fname = fullfile(pwd, 'figures', 'spin_flip_ini', subfolder, fname_str); %filename
    save_figure_pdf(h3, 10, 8, fname);
    %save_figure_eps(h3, 10, 8, fname);
end
%}
%%
K_select = [11 12 6 7];
Con_select = [20 10 25 20];
idx = zeros(5,1);

h4=figure(4);
plot(dH_avg, I_list/I_max, 'bo', 'MarkerSize', 10);
hold on
I_max = log2(N);

Con_list2 = Con_list(Con_list>0);
K_list2 = K_list(K_list>0);
for i=1:numel(K_select)
    Con_set = find(Con_list2 == Con_select(i));
    K_set = find(K_list2 == K_select(i));
    idx(i) = intersect(Con_set, K_set);
    plot(dH_avg(idx(i)), I_list(idx(i))/I_max, 'rx', 'MarkerSize', 10);
end
xlabel('$$\langle \langle dH \rangle \rangle$$');
ylabel('$$I(p_{in}, p_{out})/I_{max}$$');
set(gca, 'FontSize', 24);

qsave = 1;
dataID = 'dataset1_binomprob_summary'; % name data sets separately
if qsave
    fname_str = strrep(sprintf('%s_dH_avg_vs_information_N%d_a0_%.1f_Rcell_%.2f_flip_%d_normalized',...
        dataID, N, a0, Rcell, flip), '.', 'p');
    fname = fullfile(pwd, 'figures', 'spin_flip_ini', subfolder, fname_str); %filename
    save_figure_pdf(h4, 10, 8, fname);
    %save_figure_eps(h1, 10, 8, fname);
end
%%
Con_set = find(Con_list2 == 10);
K_set = find(K_list2 == 14);
intersect(Con_set, K_set)