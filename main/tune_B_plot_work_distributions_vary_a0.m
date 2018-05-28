%% Loads W and Delta h from from saved data to generate distributions 
% Works for protocols where a0 changes but K, Con stay constant

% Parameters 
clear variables
close all
warning off

% Parameters
gridsize = 11;
N = gridsize^2;
ntrials = 1000;
qsave = 1;
protocol = '2';
subfolder = strcat('protocol', protocol); %subfolder of ...\work_distribution to save figures in

% Load data with given tuning parameters
% (1)---Constant a0, variable Son, K---
%
a0 = 1.5;
Rcell = 0.2*a0;

% use hexagonal lattice
[dist, ~] = init_dist_hex(gridsize, gridsize);
dist_vec = a0*dist(1,:);
r = dist_vec(dist_vec>0); % exclude self influence
fN = sum(sinh(Rcell)*sum(exp(Rcell-r)./r)); % calculate signaling strength

%--- (2) constant Son, K, variable a0---
%
tsteps = 10;
a0 = linspace(1.5, 0.5, tsteps+1);

[dist, pos] = init_dist_hex(gridsize, gridsize);
dist_vec = dist(1,:);
fN = zeros(1, tsteps+1);
for i=1:tsteps+1
    r = a0(i)*dist_vec(dist_vec>0); % exclude self influence
    Rcell = 0.2*a0(i);
    fN(i) = sum(sinh(Rcell)*sum(exp(Rcell-r)./r)); % calculate signaling strength sphere
    %[out, map] = get_phenotype_map(a0, dist, Rcell, Son, K);
end

% fixed circuit parameters
Son = 5;
K = (Son+1)/2*(1+fN(1));
K_all = K*ones(1,tsteps+1);
Son_all = Son*ones(1,tsteps+1);

% list of B
B_all = (Son+1)/2*(1+fN)- K;
disp(B_all)
%

%%
% load files
name = strrep(sprintf('N%d_a0_%.1fto%.1f_K_%.1f_Son_%.f_B_%.2fto%.2f_tsteps%d_ntrials%d',...
            N, a0(1), a0(end), K, Son, B_all(1), B_all(end), tsteps, ntrials), '.', 'p');
in_file = fullfile(pwd, 'data', 'work_distribution', ...
    strcat('prot_', protocol, '_', name, '.mat'));
load(in_file);
wF = w; %forward work
qF = q;
delhF = delh;
name = strrep(sprintf('N%d_a0_%.1fto%.1f_K_%.1f_Son_%.f_B_%.2fto%.2f_tsteps%d_ntrials%d',...
            N, a0(end), a0(1), K, Son, B_all(end), B_all(1), tsteps, ntrials), '.', 'p');
in_file = fullfile(pwd, 'data', 'work_distribution', ...
    strcat('prot_', protocol, 'R_', name, '.mat'));
load(in_file); 
wR = w; % reverse work
qR = q;
delhR = delh;
%% Part I - Compare Delta h and W
h1 = figure();
hold on
line = linspace(min([wF wR])-10, max([wF wR])+10,100);
p1=plot(line, line, 'LineWidth', 1.5);
p2=scatter(wF, delhF);
p3=scatter(wR, delhR);
set(gca, 'FontSize', 24); 
set(h1, 'Units', 'inches');
set(h1, 'Position', [1 1 12 8]); 
legend([p1 p2 p3], {'W = \Delta h', 'Forward', 'Reverse'}, 'Location', 'eastoutside');
xlabel('$$W$$');
ylabel('$$\Delta h$$');
xlim([min([wF wR])-10, max([wF wR])+10]);
% save figure
if qsave
    out_file = fullfile(pwd, 'figures', 'work_distribution', subfolder,...
        strcat('W_vs_delh_prot', protocol, '_', name));
    save_figure_pdf(h1, 12, 8, out_file);
end

%% Part II - Plot P(W)
% plot separate densities
%{
% estimate densitiies
[fF, xiF] = ksdensity(wF);
%[fR, xiR] = ksdensity(wR);

% Get kernel density estimation and plot
nbins = 30;
% forward 
figure()
hold on
plot(xiF, fF, 'LineWidth', 1.5);
histogram(wF, nbins, 'Normalization', 'pdf');
% reverse
figure()
hold on
plot(xiR, fR, 'LineWidth', 1.5);
histogram(wR, nbins, 'Normalization', 'pdf');
%}

% estimate densitiies
[fF, xiF] = ksdensity(wF);
[fR, xiR] = ksdensity(wR);

% Plot forward and reverse distributions together
set(0, 'defaulttextinterpreter', 'latex');
h2 = figure();
hold on 
plot(xiF, fF, 'LineWidth', 1.5);
plot(xiR, fR, 'LineWidth', 1.5);
set(gca, 'FontSize', 24); 
set(h2, 'Units', 'inches');
set(h2, 'Position', [1 1 10 8]); 
xlabel('w');
ylabel('P(w)');
legend('Forward', 'Reverse', 'Location', 'nw');

% save figure
if qsave
    out_file = fullfile(pwd, 'figures', 'work_distribution', subfolder,...
        strcat('distribution_w_prot', protocol,'_', name));
    save_figure_pdf(h2, 10, 8, out_file);
end

%
% estimate densities in given range
%{
pts = linspace(0.4, 4.2, 100);
[fF, xiF] = ksdensity(WF, -pts);
[fR, xiR] = ksdensity(WR, pts);
ratio = fF./fR; 
ratio(ratio == Inf) = 0;
logratio = log(fF) - log(fR);
logratio(logratio == Inf) = 0;
%}

%% Part II - Plot P(Q)
% estimate densitiies
[fF, xiF] = ksdensity(qF);
[fR, xiR] = ksdensity(qR);

% Plot forward and reverse distributions together
set(0, 'defaulttextinterpreter', 'latex');
h3 = figure();
hold on 
plot(xiF, fF, 'LineWidth', 1.5);
plot(xiR, fR, 'LineWidth', 1.5);
set(gca, 'FontSize', 24); 
set(h3, 'Units', 'inches');
set(h3, 'Position', [1 1 10 8]); 
xlabel('q');
ylabel('P(q)');
legend('Forward', 'Reverse', 'Location', 'nw');

% save figure
if qsave
    out_file = fullfile(pwd, 'figures', 'work_distribution', subfolder,...
        strcat('distribution_q_prot', protocol,'_', name));
    save_figure_pdf(h3, 10, 8, out_file);
end

%% Plot ratio P_F(W)/P_R(W)
%{
h1 = figure();
hold on 
plot(pts, logratio, 'LineWidth', 1.5);
%semilogy(pts, ratio, 'LineWidth', 1.5);
%set(gca, 'FontSize', 24); 
%set(h, 'Units', 'inches');
%set(h, 'Position', [1 1 10 8]); 
xlabel('W');
ylabel('$$P_F(W)/P_R(W)$$');
%}
%% Plot Hamiltonian change P(Delta h)
[fF, xiF] = ksdensity(delhF);
[fR, xiR] = ksdensity(delhR);

h4 = figure();
hold on 
plot(xiF, fF, 'LineWidth', 1.5);
plot(xiR, fR, 'LineWidth', 1.5);
set(gca, 'FontSize', 24); 
set(h4, 'Units', 'inches');
set(h4, 'Position', [1 1 10 8]); 
xlabel('$$\Delta h$$');
ylabel('$$P(\Delta h)$$');
legend('Forward', 'Reverse', 'Location', 'ne');

% save figure
if qsave
    out_file = fullfile(pwd, 'figures', 'work_distribution', subfolder,...
        strcat('distribution_delh_prot', protocol,'_', name));
    save_figure_pdf(h4, 10, 8, out_file);
end
