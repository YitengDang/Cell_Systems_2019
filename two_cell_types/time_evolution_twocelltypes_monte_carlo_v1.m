% Time evolution of a system with visualization of the dynamics without
% noise showing the count of nearest neighbors that are ON
close all
clear all
warning off
set(0, 'defaulttextinterpreter', 'latex');
%--------------------------------------------------------------------------
% Set parameters of the system
gridsize = 15;
N = gridsize^2;
frac_type1 = 0.5;
N1 = round(frac_type1*N);
N2 = N - N1;
a0 = 1.5;
rcell = [0.2 0.2];
Rcell = rcell*a0;
Mcomm = [0 1; 1 0]; % Communication matrix, type i reacts to type j iff M_ij=1
Mcomm = logical(Mcomm);

% Circuit parameters 
Con = [15 15]; %column j = type j 
K = [2 5];
Coff = [1 1]; %default [1 1]

% Simulation parameters
p = [0.5 0.5]; %fraction of type-1 and type-2 cells which are ON
tmax = 1000;
%--------------------------------------------------------------------------
%% Calculate interaction matrix
% Distance and position
[dist, pos] = init_dist_hex(gridsize, gridsize);
[f_mat, g_mat, f_std, g_std] = calc_interaction_strengths(N, frac_type1, a0, rcell, dist, Mcomm);

if sum(sum(f_std > 10^(-2)))
    disp('Warning! Interaction strengths have large standard deviation');
end
%% Dynamics with only p1, p2
p1_all = [];
p2_all = [];
t = 0;

p1_all(end+1) = p(1);
p2_all(end+1) = p(2);

[p_out, terminate] = update_monte_carlo(N1, N2, p, Con, Coff, K, f_mat, g_mat);
while ~terminate && t<=tmax
    t = t+1;
    p1_all(end+1) = p_out(1);
    p2_all(end+1) = p_out(2);
    p = [p_out(1) p_out(2)];
    
    [p_out, terminate] = update_monte_carlo(N1, N2, p, Con, Coff, K, f_mat, g_mat);
end
%% Plot dynamics
h1 = figure(1);
hold on
plot(0:t, p1_all, 'bo-');
plot(0:t, p2_all, 'ro-');
legend({'p_1', 'p_2'});
xlabel('time');
ylabel('fraction ON cells');
ylim([0 1]);
set(gca, 'FontSize', 24);
%% Plots
%{
% hamiltonian
h2 = figure(2);
plot(0:t, h/N, 'o-');
xlabel('time');
ylabel('pseudo-energy h');
%% Plot spatial index I
h3 = figure(3);
plot(0:t, I, 'o-');
xlabel('time');
ylabel('spatial index I');
ylim([-0.2 1]);
%% Plot p1, p2
h4 = figure(4);
hold on
plot(0:t, Non1/N1, 'ro-');
plot(0:t, Non2/N2, 'bo-');
xlabel('time');
ylabel('fraction ON cells');
legend({'type 1', 'type 2'});
ylim([0 1]);
%% Plot I11, I22, I12
%{
h5 = figure(5);
hold on
plot(0:t, I_11, 'ro-');
plot(0:t, I_22, 'bo-');
plot(0:t, I_12, 'go-');
xlabel('time');
ylabel('spatial order I');
legend({'11', '22', '12'});
ylim([-0.5 1]);
%}
%% Plot theta_11, theta_22, theta_12
h6 = figure(6);
hold on
plot(0:t, theta_11, 'ro-');
plot(0:t, theta_22, 'bo-');
plot(0:t, theta_12, 'go-');
xlabel('time');
ylabel('spatial order theta');
legend({'11', '22', '12'});
ylim([-0.5 1]);
%}
%% Save trajectory
qsave = 0;
if qsave
    fname_str = strrep(...
        sprintf('N1_%d_N2_%d_p1_%.2f_p2_%.2f_a0_%.2f_K1_%d_K2_%d_Con1_%d_Con2_%d_Rcell1_%.2f_Rcell2_%.2f_t%d', ...
        N1, N2, p(1), p(2), a0, K(1), K(2), Con(1), Con(2), Rcell(1), Rcell(2), t), '.', 'p');
    i = 1;
    %fname = fullfile('K:\bn\hy\Shared\Yiteng\Multicellularity\data_twocelltypes\time_evolution', fname_str);
    fname = fullfile(pwd, 'dynamics', 'twocelltypes', ...
        strcat(fname_str,'-v',int2str(i),'.mat'));
    while exist(fname, 'file') == 2
        i=i+1;
      fname = fullfile(pwd, 'dynamics', 'twocelltypes', ...
          strcat(fname_str,'-v',int2str(i),'.mat'));
    end
    save(fname)
end