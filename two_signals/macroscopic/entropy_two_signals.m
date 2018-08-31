% Calculates the multicellular entropy for circuits with two signals
clear all
close all
clc
set(0, 'defaulttextinterpreter', 'latex');
%% Parameters
% lattice parameters
gz = 5;
N = gz^2;
a0 = 1.5;
rcell = 0.2;
Rcell = rcell*a0;

% circuit parameters 
M_int = [0 -1; 1 1];
Con = [18 16];
Coff = [1 1];
K = [0 25; 5 11];% K(i,j): sensitivity of type i to type j molecules
lambda = [1 1.2]; % diffusion length (normalize first to 1)
%hill = Inf;

% looping parameters
K12_all = [5:2:35];
K21_all = [5:2:35];
K22_all = [10];

% load dist, pos
% [dist, pos] = init_dist_hex(gz, gz);
mcsteps = 0;
[pos, dist] = initial_cells_random_markov_periodic(gz, mcsteps, rcell);

% Calculate interaction strength
dist_vec = a0*dist(1,:);
r = dist_vec(dist_vec>0); % exclude self influence
fN = zeros(1, 2);
gN = zeros(1, 2);
for i=1:2
    fN(i) = sum(sinh(Rcell) * sum(exp((Rcell-r)./lambda(i)).*(lambda(i)./r)) ); % calculate signaling strength
    gN(i) = sum(sinh(Rcell)^2 * sum(exp(2*(Rcell-r)./lambda(i)).*(lambda(i)./r).^2 ) ); % calculate noise variance strength
end

%% Tester code
%{
% fix n, calculate S_{E,n}
n = [0 5; 1 19];
pij = n/N;

W = transition_prob_states_two_signals_analytical_calc(M_int, Con, Coff, K, fN, gN, pij);
iniON = reshape(round(pij.*N)',4,1); % iniON: (0,0), (0,1), (1,0), (1,1)
Peq = prod(diag(W).^iniON); % diag(W): W(0,0), W(0,1), W(1,0), W(1,1)
disp(Peq)

% multinomial coefficient
mult_coeff = prod([1:N]./[1:n(1) 1:n(2) 1:n(3) 1:n(4)]); 
% mult_coeff_2 = nchoosek(n(1)+n(2), n(2))*nchoosek(n(1)+n(2)+n(3),
% n(3))*nchoosek(n(1)+n(2)+n(3)+n(4), n(4)); % alternative

Omega_En = mult_coeff*Peq;

% sum over all n
[Omega_E, S, frac_E] = entropy_two_signals_calc(N, M_int, Con, Coff, K, fN, gN);
%}
%% Loop over parameter values
S_all = zeros(numel(K12_all), numel(K21_all), numel(K22_all));
for i1=1:numel(K12_all)
    for i2=1:numel(K21_all)
        for i3=1:numel(K22_all)
            fprintf('Loop indices %d %d %d \n', i1, i2, i3);
            this_K = K;
            this_K(1,2) = K12_all(i1);
            this_K(2,1) = K21_all(i2);
            this_K(2,2) = K22_all(i3);
            [~, S, ~] = entropy_two_signals_calc(N, M_int, Con, Coff, this_K, fN, gN);
            %Omega_E_all(i1,i2,i3) = Omega_E;
            S_all(i1,i2,i3) = S;
        end
    end
end

Smax = N*log(2);
frac_E_all = exp(S_all - Smax);
%% Load results
%{
folder = 'H:\My Documents\Multicellular automaton\figures\two_signals\entropy';
fname_str = 'Entropy_M_int_0_-1_1_1_N25_a0_1p50_rcell_0p20_K_sweep_all_Con_18_16';
load(fullfile(folder, fname_str));
%}
%% Plot results
K22_idx = 1; % fix K22
this_K = K; this_K(2,2) = K22_all(K22_idx);

% Plot S
h=figure;
imagesc(K12_all, K21_all, S_all(:,:,K22_idx)' );
hold on
%for i=1:3
%    line([0.5, 4.5], [i+0.5, i+0.5], 'Color', 'r', 'LineWidth', 1);
%    line([i+0.5, i+0.5], [0.5, 4.5], 'Color', 'r', 'LineWidth', 1);
%end
c = colorbar;
set(gca, 'Ydir', 'normal', 'FontSize', 20)
xlabel('$$K^{(12)}$$', 'FontSize', 24)
ylabel('$$K^{(21)}$$', 'FontSize', 24)
c.Label.String = 'S';
%set(gca, 'XTick', K12_all, 'YTick', K21_all);
set(gca, 'XTick', K12_all(1):5:K12_all(end), 'YTick', K21_all(1):5:K21_all(end));
set(h, 'Position', [100 100 640 520]);
%caxis([0 1]);

%% Plot frac equilibrium
h2 = figure;
imagesc(K12_all, K21_all, frac_E_all(:,:,K22_idx)' );
hold on
c = colorbar;
%caxis([0 1]);
title('Fraction equilibrium states');
c.Label.String = 'Fraction';
set(gca, 'Ydir', 'normal', 'FontSize', 20)
xlabel('$$K^{(12)}$$', 'FontSize', 24)
ylabel('$$K^{(21)}$$', 'FontSize', 24)
set(gca, 'XTick', K12_all(1):5:K12_all(end), 'YTick', K21_all(1):5:K21_all(end));
set(h2, 'Position', [100 100 640 520]);

%% Save figures
% data strings
M_int_s = sprintf('%d_%d_%d_%d', M_int(1,1), M_int(1,2), M_int(2,1), M_int(2,2));
a0_s = sprintf('%.2f', a0);
R_s = sprintf('%.2f', rcell);
K_s = sprintf('%d_x_y_%d', this_K(1,1), this_K(2,2));
Con_s = sprintf('%d_%d', Con(1), Con(2));

folder = 'H:\My Documents\Multicellular automaton\figures\two_signals\entropy';
fname_str = strrep(sprintf('Entropy_M_int_%s_N%d_a0_%s_rcell_%s_K_%s_Con_%s',...
    M_int_s, N, a0_s, R_s, K_s, Con_s), '.', 'p');
fname_str_2 = strrep(sprintf('Entropy_M_int_%s_N%d_a0_%s_rcell_%s_K_%s_Con_%s',...
    M_int_s, N, a0_s, R_s, 'sweep_all', Con_s), '.', 'p');

% Save data
qsave = 0;
if qsave
save(fullfile(folder, 'data', fname_str_2), 'gz', 'N', 'a0', 'rcell', ...
    'Rcell', 'M_int', 'Con', 'Coff', 'K', 'lambda', 'K12_all', 'K21_all', 'K22_all',...
    'S_all', 'frac_E_all');
end

% Save figure
qsave = 0;
save_figure(h, 10, 8, fullfile(folder, fname_str), '.pdf', qsave);
save_figure(h2, 10, 8, fullfile(folder, strcat(fname_str, '_frac_eq')), '.pdf', qsave);
