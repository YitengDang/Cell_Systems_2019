% Obtains the multicellular entropy for circuits with two signals
% from simulations
clear all
close all
clc
set(0, 'defaulttextinterpreter', 'latex');

%% Likely to work only for small enough N
%{
gz = [5:25];
N = gz.^2;
k = 4;
binom = zeros(numel(N), 1);

for i=1:numel(N)
    binom(i) = nchoosek(N(i)+k-1, k);
end

figure;
plot(gz, binom, 'o-');
xlabel('grid size');
ylabel('# points');
title('# points on standard 4-simplex');
set(gca, 'YScale', 'log');
%}
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
K = [0 25; 5 10];% K(i,j): sensitivity of type i to type j molecules
lambda = [1 1.2]; % diffusion length (normalize first to 1)
hill = Inf;
noise = 0;

% looping parameters
K12_all = [5:5:35];
K21_all = [5:5:35];
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
n = [0 0; 24 1];
%pij = n/N;

% multinomial coefficient
mult_coeff = prod([1:N]./[1:n(1) 1:n(2) 1:n(3) 1:n(4)]); 
% mult_coeff_2 = nchoosek(n(1)+n(2), n(2))*nchoosek(n(1)+n(2)+n(3),
% n(3))*nchoosek(n(1)+n(2)+n(3)+n(4), n(4)); % alternative

% Simulate P_eq
max_sim = 10^4;
Q = 0; % number of equilibrium states
nsim = min(max_sim, mult_coeff);
fprintf('Simulations = %d \n', nsim);
for i=1:nsim
    disp(i);
    % initiate random config
        % subdivide indices for different cell states
    idx = cell(3,1);
    idx_not00 = randperm(N, N-n(1,1)); % cell indices of all cells not (0,0)
    part_idx = [1 n(1,2) n(1,2)+n(2,1) n(1,2)+n(2,1)+n(2,2)]; % partition indices
    idx{1} = idx_not00(part_idx(1):part_idx(2));
    idx{2} = idx_not00(part_idx(2)+1:part_idx(3));
    idx{3} = idx_not00(part_idx(3)+1:part_idx(4));

    %{
    %Check that sets have the right number of elems
    check = (numel(idx{1})==iniON(1,2));
    check = check*(numel(idx{2})==iniON(2,1));
    check = check*(numel(idx{3})==iniON(2,2));
    % Check that sets have disjoint elements
    check = check*(numel(union(union(idx{1}, idx{2}), idx{3}))==N-iniON(1,1)); 
    fprintf('Check passed? %d \n', check);
    %}

    % turn genes ON
    cells = zeros(N, 2);
    cells(idx{1}, 2) = 1;
    cells(idx{2}, 1) = 1;
    cells(idx{3}, 1) = 1;
    cells(idx{3}, 2) = 1;
    
    % check if stable
    [~, changed] = update_cells_two_signals_multiply_finite_Hill(...
        cells, dist, M_int, a0, Rcell, Con, Coff, K, lambda, hill, noise);
    if ~changed
        Q=Q+1;
    end
end
Peq = Q/nsim;

Omega_En = mult_coeff*Peq;

% sum over all n
%[Omega_E, S, frac_E] = entropy_two_signals_calc_sim(N, M_int, Con, Coff, K, fN, gN);
%}
%%
%{
this_K = [0 5; 5 10];
entropy_two_signals_sim_make_dir(gz, a0, rcell, M_int, Con, this_K, lambda, hill, noise)
[~, S, ~] = entropy_two_signals_calc_sim(N, a0, rcell, lambda,...
    hill, noise, M_int, Con, Coff, this_K, dist);
%}
%% Loop over parameter values
%parent_folder = 'H:\My Documents\Multicellular automaton\data\two_signals\entropy'; % for saving
parent_folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\entropy';

S_all = zeros(numel(K12_all), numel(K21_all), numel(K22_all));
for i1=1:numel(K12_all)
    for i2=1:numel(K21_all)
        for i3=1:numel(K22_all)
            fprintf('Loop indices %d %d %d \n', i1, i2, i3);
            this_K = K;
            this_K(1,2) = K12_all(i1);
            this_K(2,1) = K21_all(i2);
            this_K(2,2) = K22_all(i3);
            
            entropy_two_signals_sim_make_dir(parent_folder, ...
                gz, a0, rcell, M_int, Con, this_K, lambda, hill, noise)

            [~, S, ~] = entropy_two_signals_calc_sim(N, a0, rcell, lambda,...
                hill, noise, M_int, Con, Coff, this_K, dist, parent_folder);
            %Omega_E_all(i1,i2,i3) = Omega_E;
            S_all(i1,i2,i3) = S;
        end
    end
end
%}
Smax = N*log(2);
frac_E_all = exp(S_all - Smax);
%}

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
fname_str = strrep(sprintf(...
    'Sim_entropy_M_int_%s_N%d_a0_%s_rcell_%s_K_%s_Con_%s_hill%.1f_noise%.1f',...
    M_int_s, N, a0_s, R_s, K_s, Con_s, hill, noise), '.', 'p');
fname_str_2 = strrep(sprintf(...
    'Sim_entropy_M_int_%s_N%d_a0_%s_rcell_%s_K_%s_Con_%s_hill%.1f_noise%.1f',...
    M_int_s, N, a0_s, R_s, 'sweep_all', Con_s, hill, noise), '.', 'p');

% Save data
qsave = 1;
if qsave
save(fullfile(folder, 'data', fname_str_2), 'gz', 'N', 'a0', 'rcell', ...
    'Rcell', 'M_int', 'Con', 'Coff', 'K', 'lambda', 'K12_all', 'K21_all', 'K22_all',...
    'S_all', 'frac_E_all');
end

% Save figure
qsave = 1;
save_figure(h, 10, 8, fullfile(folder, fname_str), '.pdf', qsave);
save_figure(h2, 10, 8, fullfile(folder, strcat(fname_str, '_frac_eq')), '.pdf', qsave);