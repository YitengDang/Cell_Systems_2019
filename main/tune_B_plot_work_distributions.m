%% Loads W and Delta h from from saved data to generate distributions 
% Works for protocols where a0 does not change

% Parameters 
clear variables
close all
warning off
clc

%% Specify parameters
gridsize = 11;
N = gridsize^2;
ntrials = [1000 1000]; % (1) forward, (2) reverse
qsave = 0;
protocol = 4;
protid = num2str(protocol);
subfolder = strcat('protocol', protid);

% (1)---Constant a0, variable Son, K---
%
a0 = 1.5;
Rcell = 0.2*a0;

% use hexagonal lattice
[dist, pos] = init_dist_hex(gridsize, gridsize);
dist_vec = a0*dist(1,:);
r = dist_vec(dist_vec>0); % exclude self influence
fN = sum(sinh(Rcell)*sum(exp(Rcell-r)./r)); % calculate signaling strength

% protocols 1, 1v2(3)
%{
tsteps = 10; % complete procedure in tsteps steps
Son_0 = 5; 
K_ini = (Son_0+1)/2*(1+fN);
if protocol == 3
    load('tuneB_protocol3_K4p5to4p5_Son5to11p5_B0to5.mat');
else
    Son = linspace(5, 15, tsteps+1);
end
K = linspace(K_ini, K_ini, tsteps+1);
%}
% protocol 4 & 6
%
tsteps = 10;
K_ini = 12;
Con_ini = 2*K_ini/(1+fN) - 1;
K = linspace(K_ini,2,tsteps+1); %protocol 4
%K = linspace(K_ini,4,tsteps+1); %protocol 6
Con = linspace(Con_ini, Con_ini, tsteps+1);
%}

B_all = (Con+1)/2*(1+fN) - K;
disp('B='); disp(B_all)
%% Load data
if find(protocol == [1 3 4 5 6], 1)
    name = strrep(sprintf('N%d_a0_%.1f_K_%.1fto%.1f_Son_%.fto%.f_B_%.2fto%.2f_tsteps%d_ntrials%d', ...
                N, a0, K(1), K(end), Con(1), Con(end), B_all(1), B_all(end),...
                tsteps, ntrials(1)), '.', 'p');     
elseif protocol == 2
    name = strrep(sprintf('N%d_a0_%.1fto%.1f_K_%.1f_Son_%.f_B_%.2fto%.2f_tsteps%d_ntrials%d', ...
            N, a0(1), a0(end), K, Con, B_all(1), B_all(end),...
            tsteps, ntrials(1)), '.', 'p');
else 
    disp('ERROR');
end

% Forward trajectories
in_file = fullfile(pwd, 'data', 'work_distribution', ...
    strcat('prot_', protid,'_', name, '.mat'));
load(in_file);
wF = w; %forward work
qF = q;
delhF = delh;
%%
% Reverse trajectories
if find(protocol == [1 3 5], 1)
    name = strrep(sprintf('N%d_a0_%.1f_K_%.1fto%.1f_Son_%.fto%.f_B_%.2fto%.2f_tsteps%d_ntrials%d', ...
                N, a0, K(1), K(end), Con(end), Con(1), B_all(end), B_all(1), tsteps, ntrials(2)), '.', 'p'); 
elseif find(protocol == [4 6], 1)
        name = strrep(sprintf('N%d_a0_%.1f_K_%.1fto%.1f_Son_%.fto%.f_B_%.2fto%.2f_tsteps%d_ntrials%d', ...
                N, a0, K(end), K(1), Con(1), Con(end), B_all(end), B_all(1), tsteps, ntrials(2)), '.', 'p'); 
elseif protocol == 2
    name = strrep(sprintf('N%d_a0_%.1fto%.1f_K_%.1f_Son_%.f_B_%.2fto%.2f_tsteps%d_ntrials%d', ...
            N, a0(1), a0(end), K, Con, B_all(1), B_all(end),...
            tsteps, ntrials(2)), '.', 'p');
else
    disp('Error!');
end

in_file = fullfile(pwd, 'data', 'work_distribution', ...
    strcat('prot_', protid, 'R_', name, '.mat'));
load(in_file);
wR = w; %reverse work
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
        strcat('W_vs_delh_prot', protid, '_', name));
    save_figure_pdf(h1, 12, 8, out_file);
end

%% II - Plot P(W)
% estimate densitiies
[fF, xiF] = ksdensity(wF);
%[fR, xiR] = ksdensity(wR);

% Plot forward and reverse distributions together
set(0, 'defaulttextinterpreter', 'latex');
h2 = figure();
hold on 
plot(xiF, fF, 'LineWidth', 1.5);
%plot(xiR, fR, 'LineWidth', 1.5);
set(gca, 'FontSize', 24); 
set(h2, 'Units', 'inches');
set(h2, 'Position', [1 1 10 8]); 
xlabel('w');
ylabel('P(w)');
%legend('Forward', 'Reverse', 'Location', 'nw');

% save figure without means
if qsave
    out_file = fullfile(pwd, 'figures', 'work_distribution', subfolder,...
        strcat('distribution_w_prot', protid,'_', name));
    save_figure_pdf(h2, 10, 8, out_file);
end

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
%
% Compare with <W> from distribution and plot
dW = xiF(2)-xiF(1); %step size 
W_av = sum(xiF.*fF)*dW;
h2;
plot([W_av W_av], [0 round(max(fF), 1)]);

%% Compare with <W> from transition matrix
% Load data
% protocol 4
fileid = strrep(sprintf('wpt_N%d_Con_%.1f_K%.1fto%.1f_a0_%.1f_Rcell%.1f',...
    N, Con_ini, K(1), K(end), a0, Rcell), '.', 'p');
in_file = fullfile(pwd, 'data', 'distribution_pI_from_tmatrix', fileid);
load(in_file, 'W');

% Plot <W> in same figure as P(W)
h2;
plot([W W], [0 round(max(fF), 1)]);
legend('Forward', 'Distr. mean', 'T-matrix mean', 'Location', 'nw');

% save figure
if qsave
    out_file = fullfile(pwd, 'figures', 'work_distribution', subfolder,...
        strcat('distribution_w_avgw_prot', protid,'_', name));
    save_figure_pdf(h2, 10, 8, out_file);
end

%% III - Plot P(Q)
% estimate densitiies
[fF, xiF] = ksdensity(qF);
%[fR, xiR] = ksdensity(qR);

% Plot forward and reverse distributions together
set(0, 'defaulttextinterpreter', 'latex');
h3 = figure();
hold on 
plot(xiF, fF, 'LineWidth', 1.5);
%plot(xiR, fR, 'LineWidth', 1.5);
set(gca, 'FontSize', 24); 
set(h3, 'Units', 'inches');
set(h3, 'Position', [1 1 10 8]); 
xlabel('q');
ylabel('P(q)');
legend('Forward', 'Reverse', 'Location', 'nw');

% save figure
if qsave
    out_file = fullfile(pwd, 'figures', 'work_distribution', subfolder,...
        strcat('distribution_q_prot', protid,'_', name));
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
%% IV - Plot Hamiltonian change P(Delta h)
[fF, xiF] = ksdensity(delhF);
%[fR, xiR] = ksdensity(delhR);

h4 = figure();
hold on 
plot(xiF, fF, 'LineWidth', 1.5);
%plot(xiR, fR, 'LineWidth', 1.5);
set(gca, 'FontSize', 24); 
set(h4, 'Units', 'inches');
set(h4, 'Position', [1 1 10 8]); 
xlabel('$$\Delta h$$');
ylabel('$$P(\Delta h)$$');
legend('Forward', 'Reverse', 'Location', 'ne');

% save figure
if qsave
    out_file = fullfile(pwd, 'figures', 'work_distribution', subfolder,...
        strcat('distribution_delh_prot', protid,'_', name));
    save_figure_pdf(h4, 10, 8, out_file);
end