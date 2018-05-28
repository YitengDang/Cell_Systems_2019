% Calculates the number of accessible states from a given value of the
% Hamiltonian
clear all
close all
warning off
clc
%% Parameters
gridsize = 11;
N = gridsize^2;
a0 = 0.5;
Rcell = 0.2*a0;
K = 16;
Con = 8;

% Calculate fN
[pos,ex,ey] = init_cellpos_hex(gridsize,gridsize);
dist = dist_mat(pos,gridsize,gridsize,ex,ey);
dist_vec = a0*dist(1,:);
r = dist_vec(dist_vec>0); % exclude self influence
fN = sum(sinh(Rcell)*sum(exp(Rcell-r)./r)); % calculate signaling strength

%% Number of accessible states from DOS
% load DOS data
% if error, check other parameters in code
path = fullfile(pwd, 'data', 'dos', 'probI_fixedp','2017-06-01'); %directory to look up data
[log_omegapI, pv, Iv] = dos_Wang_Landau_omegapI_calc(N, a0, path);
%fname_str = sprintf('OmegapI_WL_norm_N121_a0_1p50_f0exp1_ffin_e10e-4_pflat0p8_max_39_bins.mat');
%fname = fullfile(pwd, 'data', 'dos', 'OmegapI_probI_fixedp', fname_str); %directory to look up data
%load(fname, 'omegapI_all', 'n', 'Icenters');

% Plot DOS
set(0, 'defaulttextinterpreter', 'latex');
h1=figure(1);
colormap('hot');
%imagesc(pv, Iv, log_OmegapI);
if mod(N,2)
    %pv = [n/N (N+1-n(1:end-1))/N];
    %Iv = Icenters;
    %log_omegapI = log([omegapI_all  fliplr(omegapI_all(:,1:end-1))]);
    imagesc(pv, Iv, log_omegapI);
else
    %pv = [n/N (N+1-n)/N];
    %Iv = Icenters;
    %log_omegapI = log([omegapI_all  fliplr(omegapI_all(:,1:end))]);
	imagesc(pv, Iv, log_omegapI);
end
c=colorbar;
c.Label.String = 'log(\Omega(p,I))';
set(gca,'FontSize', 24);
set(gca,'YDir','normal');
xticks(0.1:0.1:0.9);
%yticks(edges(1:3:end));
xlabel('p');
ylabel('I');
%}
%% Plot contour map of h
h = @(p, I) -0.5*(Con-1)*(1 + 4*fN.*p.*(1-p).*I + fN*(2*p-1).^2)...
    -(2*p-1).*(0.5*(Con+1)*(1+fN) - K);
%pv = (0:N)/N;
%Iv = -0.05:0.01:1;
%pv = n/N;
pv = sort(pv);
[pmesh, Imesh] = meshgrid(pv, Iv);
%[pmesh, Imesh] = meshgrid(piv, Iv);
hmesh = h(pmesh, Imesh);

% Calculate maximum h
hmax = max(max(hmesh));
hmin = min(min(hmesh));
[i1, i2] = find(hmesh==hmax);
[j1, j2] = find(hmesh==hmin);

h2=figure(2);
hold on
contourf(pv, Iv, hmesh, 'LineStyle', 'none');
scatter(pv(i2), Iv(i1), 'r*'); %maximum h
scatter(pv(j2), Iv(j1), 'b*'); %maximum h
colormap('summer')
c = colorbar;
ylim(c, [floor(hmin) ceil(hmax)])
xlabel('p');
ylabel('I');
ylabel(c, 'h');
set(gca, 'FontSize', 24);

qsave=0;
if qsave
    fname_str = strrep(sprintf('Hamiltonian_map_N%d_a0_%.2f_K%d_Con%d', N,...
        a0, K, Con), '.', 'p');
    fname = fullfile(pwd, 'rebuttal', 'accessible_states', fname_str);
    save_figure_pdf(h2, 10, 8, fname);
end
%% Accessible region
%{
% Specify h to consider
this_h = -4;
h_allowed = hmesh<this_h;

h3=figure(3);
hold on
imagesc(pv, Iv, h_allowed);
set(gca, 'YDir', 'normal');
colormap('gray')
%c = colorbar;

% plot Imax
fa0 = sinh(Rcell)*exp(Rcell-a0)/a0;
Imax = @(p) (6 - 4./sqrt(p*N) - 6./(p*N) - 6*p)*fa0./(1-p)/fN;
Imax2 = @(p) (6*p - 4./sqrt((1-p)*N) - 6./((1-p)*N))*fa0./p/fN;
plot(pv, max(Imax(pv), Imax2(pv)), '--b', 'Linewidth', 1.5)
xlim([0 1]);
ylim([-0.05 ceil(10*max(max(Imax(pv), Imax2(pv))))/10]);
set(gca, 'FontSize', 24);
%}
%% Plot number of accessible states (Omega_A) against h
hlist = -4:-1:-14;
astates = zeros(numel(hlist), 1);
count = 0;
% Set invalid omegapI values to 0
idx = find(log_omegapI==Inf);
log_omegapI(idx) = 0;
for h = hlist
    count = count+1;
    h_allowed = hmesh<h;
    astates(count, 1) = log(sum(sum(h_allowed.*exp(log_omegapI))));
    %astates(count, 1) = sum(sum(h_allowed.*omegapI_all));
end

h4=figure(4);
plot(hlist, astates, 'o-', 'LineWidth', 2);
xlabel('$$h$$');
ylabel('$$\log(\Omega_A)$$');
set(gca, 'FontSize', 24);
%xlim([-9 -4]);

qsave=0;
if qsave
    fname_str = strrep(sprintf('OmegaA_vs_h_N%d_a0_%.2f_K%d_Con%d', N,...
        a0, K, Con), '.', 'p');
    fname = fullfile(pwd, 'rebuttal', 'accessible_states', fname_str);
    save_figure_pdf(h4, 10, 8, fname);
end
%% Plot dynamics of a trajectory in terms of number of accessible states
%
% Load sample trajectory
dir = 'H:\My Documents\Multicellular automaton\data\dynamics\no_noise'; 
fname = fullfile(dir, 'N121_n24_neq_121_a05_K_10_Son_21_t_2-v4.mat');
load(fname, 'mom', 't', 'p');

% Calculate # accessible states along trajectory
astates_traj = zeros(t+1, 1);
h_traj = zeros(t+1, 1);
for i=1:t+1
    h_traj(i) = mom(i)/N;
    h_allowed = hmesh<h_traj(i);
    astates_traj(i, 1) = log(sum(sum(h_allowed.*exp(log_omegapI))));
end

% Plot number of accessible states
h5=figure(5);
hold on
yyaxis left
%plot(h_traj, astates_traj, 'o-'); % against h
plot(0:t, astates_traj, 'o-', 'LineWidth', 2); % against time
xlabel('$$t$$');
ylabel('$$\log(\Omega_A)$$');

yyaxis right
plot(0:t, mom/N, '-o', 'LineWidth', 2);
ylabel('$$h$$');
set(gca, 'FontSize', 24);
%xlim([-9 -4]);

% Save trajectory
qsave=1;
if qsave
    fname_str = strrep(sprintf('OmegaA_vs_t_N%d_a0_%.2f_K%d_Con%d_p%.1f_t%d', N,...
        a0, K, Con, p, t), '.', 'p');
    fname = fullfile(pwd, 'rebuttal', 'accessible_states', fname_str);
    save_figure_pdf(h5, 10, 8, fname);
end
%% Save all data
qsave=0;
if qsave
fname_str = strrep(sprintf('Accessible_states_N%d_a0_%.2f_K%d_Con%d', N, a0, K, Con), '.', 'p');
fname = fullfile(pwd, 'data', 'accessible_states', strcat(fname_str, '.mat'));
save(fname);
end
%}