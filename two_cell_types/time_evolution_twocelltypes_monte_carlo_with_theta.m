% Simulates the Monte Carlo dynamics of the system with two cell types with
% spatial order
close all
clear variables
%warning off
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
Mcomm = [1 1; 1 1]; % Communication matrix, type i reacts to type j iff M_ij=1
Mcomm = logical(Mcomm);

% Circuit parameters
Con = [8 8]; %column j = type j 
K = [10 10];
Coff = [1 1]; %default [1 1]

% Simulation parameters
p_in = [0.5; 0.5]; %fraction of type-1 and type-2 cells which are ON
% p_in must be column vector!
tmax = 1000;
%--------------------------------------------------------------------------
%% Calculate interaction matrix
% Distance and position
[dist, pos] = init_dist_hex(gridsize, gridsize);
[f_mat, g_mat, f_std, g_std] = calc_interaction_strengths(N, frac_type1, a0, rcell, dist, Mcomm);
%f_mat = [f11 f12; f21 f22];
%g_mat = [g11 g12; g21 g22];
if sum(sum(f_std > 10^(-2)))
    disp('Warning! Interaction strengths have large standard deviation');
end

%% Dynamics with only p1, p2
p1_all = [];
p2_all = [];
theta11_all = [];
theta12_all = [];
theta22_all = [];
pe1 = [];
pe2 = [];

t = 0;
theta_in = diag(2*p_in-1)*f_mat*diag(2*p_in-1);
%%
p1_all(1) = p_in(1);
p2_all(1) = p_in(2);
theta11_all(1) = theta_in(1,1);
theta12_all(1) = theta_in(1,2);
theta22_all(1) = theta_in(2,2);
[p_out, theta_out, pe, cont] = ...
    update_monte_carlo_with_theta(N1, N2, p_in, theta_in, Con, Coff, K, f_mat, g_mat);
pe1(1) = pe(1);
pe2(1) = pe(2);

while cont && t<tmax
    t = t+1;
    p1_all(end+1) = p_out(1);
    p2_all(end+1) = p_out(2);
    theta11_all(end+1) = theta_out(1,1);
    theta12_all(end+1) = theta_out(1,2);
    theta22_all(end+1) = theta_out(2,2);
    
    p_in = p_out;
    theta_in = theta_out;
    [p_out, theta_out, pe, cont] = ...
        update_monte_carlo_with_theta(N1, N2, p_in, theta_in, Con, Coff, K, f_mat, g_mat);
    pe1(end+1) = pe(1);
    pe2(end+1) = pe(2);
end
%% Plot dynamics
% Plot p1, p2
h1 = figure(1);
hold on
plot(0:t, p1_all(1:t+1), 'bo-');
plot(0:t, p2_all(1:t+1), 'ro-');
legend({'p_1', 'p_2'});
xlabel('time');
ylabel('fraction ON cells');
ylim([0 1]);
set(gca, 'FontSize', 24);
set(h1, 'Units', 'Inches', 'Position', [1 1 10 8]);
%%
% Plot theta
h2 = figure(2);
plot(0:t, theta11_all, 'b.-'); % termination probability
hold on
plot(0:t, theta12_all, 'r.-'); % termination probability
plot(0:t, theta22_all, 'g.-'); % termination probability
legend({'theta_{11}', 'theta_{12}', 'theta_{22}'}, 'Location', 'nw');
xlabel('time');
ylabel('spatial order');
%xlim([0 t]);
ylim([-1 1]);
set(gca, 'FontSize', 24);
set(h2, 'Units', 'Inches', 'Position', [1 1 10 8]);
%%
% Plot equilibrium probabilities
h3 = figure(3);
hold on
plot(0:t, pe1, 'b.-'); % termination probability
plot(0:t, pe2, 'r.-'); % termination probability
legend({'p_{eq,1}', 'p_{eq,2}'}, 'Location', 'nw');
xlabel('time');
ylabel('eq. probability');
%xlim([0 t]);
ylim([0 1]);
set(gca, 'FontSize', 24);
set(h3, 'Units', 'Inches', 'Position', [1 1 10 8]);
%% Save trajectory
qsave = 0;
if qsave
    fname_str = strrep(...
        sprintf('MC_with_theta_N1_%d_N2_%d_p1_%.2f_p2_%.2f_a0_%.2f_K1_%d_K2_%d_Con1_%d_Con2_%d_Rcell1_%.2f_Rcell2_%.2f_t%d', ...
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