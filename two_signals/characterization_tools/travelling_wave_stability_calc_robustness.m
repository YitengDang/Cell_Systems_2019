clear all
close all
set(0, 'defaulttextinterpreter', 'tex')

%% Setup
% parameters
network = 15;
gz = 15;
N = gz^2;
a0 = 1.5;
rcell = 0.2;
lambda = [1 1.2];

% get pos, dist
mcsteps = 0;
[~, dist] = initial_cells_random_markov_periodic(gz, mcsteps, rcell);
[fN, fnn, ~] = calc_fN(a0, rcell, gz, dist, lambda);

% mean-field 
pw = [1/8 1/8];
Rcell = rcell*a0;
f_MF = (4*pi/a0^2/sqrt(3)).*exp(Rcell./lambda).*sinh(Rcell./lambda).*exp(-3*a0./(2*lambda));

L = 1000;
N_inf = 0;

% save figure folder
save_folder = 'H:\My Documents\Multicellular automaton\figures\trav_wave_stability\robustness_analytical';

% Test case
%{
gz = 15;
[~, dist] = initial_cells_random_markov_periodic(gz, mcsteps, rcell);
a0 = 0.1;
[fN, fnn, ~] = calc_fN(a0, rcell, gz, dist, lambda);

[Q_value] = calc_Q_value(L, fN, fnn, pw, N_inf, f_MF);
%}
% Plot region boundaries (check)
%{
showlines = 1;
h = plot_boundaries_K_Con(fN, fnn, pw, L, showlines);

qsave = 0;
fname_str = 'plot_boundaries_test_case_1';
fname = fullfile(save_folder, fname_str);
save_figure(h, 10, 8, fname, '.pdf', qsave);
%}
%% ----------Study effect of varying parameters----------------------------
% Vary gz
gz_all = 15; %5:50; %[5 15 50]; %[5 15 50]; %[5 15 50]; %5:5:50;
% Vary a0
a0_all = 1.5; %0.1:0.1:5; %[0.5 1.5 5]; %0.1:0.01:5; %0.1:0.1:5;
% Vary lambda1, lambda2 
%lambda1_all = 1;
%lambda2_all = 1.2; %4; %:5;

% Vary l1, l2 (l = lambda/a0)
l1_all = 0.1:0.1:3;
l2_all = 0.1:0.1:4;

% results to store
% --------- choose between lambda and l -----------------------------------
%Q_value_all = zeros(numel(gz_all), numel(a0_all), numel(lambda1_all), numel(lambda2_all));
Q_value_all = zeros(numel(gz_all), numel(a0_all), numel(l1_all), numel(l2_all));

for j=1:numel(gz_all)
    this_gz = gz_all(j);
    % Calculate distances
    [~, this_dist] = initial_cells_random_markov_periodic(this_gz, mcsteps, rcell);
    
    for i=1:numel(a0_all)
        this_a0 = a0_all(i);
        %fprintf('gz = %d, a0 = %.1f \n', this_gz, this_a0);
        for i1=1:numel(l1_all) %numel(lambda1_all) %
            for i2=1:numel(l2_all) %numel(lambda2_all) %
                % -------Option 1: use lambda------------------------------
                %{
                this_lambda = [lambda1_all(i1) lambda2_all(i2)];
                fprintf('gz = %d, a0 = %.1f, lambda = [%.1f %.1f] \n',...
                    this_gz, this_a0, this_lambda(1), this_lambda(2));
                %}
                % -------Option 2: use l-----------------------------------
                %
                this_l = [l1_all(i1) l2_all(i2)];
                this_lambda = this_l*this_a0;
                fprintf('gz = %d, a0 = %.1f, l = [%.1f %.1f] \n',...
                    this_gz, this_a0, this_l(1), this_l(2));
                %}                
                % -----------Calculate interaction strengths-----------------------
                [fN, fnn, ~] = calc_fN(this_a0, rcell, this_gz, this_dist, this_lambda);                
                %{
                Rcell = rcell*this_a0;
                f_MF = (4*pi/this_a0^2/sqrt(3)).*exp(Rcell./this_lambda).*sinh(Rcell./this_lambda).*exp(-3*this_a0./(2*this_lambda));
                %}
                %------------------------------------------------------------------
                % (1) Fixed phase space region [1,L]x[1,L]
                L = 1000;
                N_inf = 0;
                Q_value = calc_Q_value(L, fN, fnn, pw, N_inf, f_MF);
                %------------------------------
                % Store results
                Q_value_all(j, i, i1, i2) = Q_value;
            end
        end
    end
end

% Save data 
% l1, l2
%
fname_str = strrep(sprintf('Q_values_gz_%d_to_%d_a0_%.1f_to_%.1f_l1_%.1f_to_%.1f_l2_%.1f_to_%.1f',...
    gz_all(1), gz_all(end), a0_all(1), a0_all(end), l1_all(1), l1_all(end), l2_all(1), l2_all(end)), '.', 'p');
fname = fullfile( save_folder, 'data', fname_str );
save( fname, 'gz_all', 'a0_all', 'l1_all', 'l2_all', 'Q_value_all' );
%}
%{
% lambda1, lambda2
fname_str = strrep(sprintf('Q_values_gz_%d_to_%d_a0_%.1f_to_%.1f_l1_%.1f_to_%.1f_l2_%.1f_to_%.1f',...
    gz_all(1), gz_all(end), a0_all(1), a0_all(end), lambda1_all(1), lambda1_all(end), lambda2_all(1), lambda2_all(end)), '.', 'p');
fname = fullfile( save_folder, 'data', fname_str );
save( fname, 'gz_all', 'a0_all', 'lambda1_all', 'lambda2_all', 'Q_value_all' );
%}
%% Load saved data
fname_str = strrep(sprintf('Q_values_gz_%d_to_%d_a0_%.1f_to_%.1f_l1_%.1f_to_%.1f_l2_%.1f_to_%.1f',...
    gz_all(1), gz_all(end), a0_all(1), a0_all(end), l1_all(1), l1_all(end), l2_all(1), l2_all(end)), '.', 'p');
fname = fullfile( save_folder, 'data', fname_str );
load( fname, 'gz_all', 'a0_all', 'l1_all', 'l2_all', 'Q_value_all' );

%% Plot and save examples of phase plots
%{
pw = [1/8 1/8]; % Fix pw
L = 1000;
this_gz = 15;
this_a0 = 0.2;
this_lambda = lambda;

[~, this_dist] = initial_cells_random_markov_periodic(this_gz, mcsteps, rcell);
[fN, fnn, ~] = calc_fN(this_a0, rcell, this_gz, this_dist, this_lambda);
%}
this_gz = 15;
this_a0 = 1;
this_l = [2.2 100]; %[2.2 0.7];
this_lambda = [this_l(1) this_l(2)]*a0;
[~, this_dist] = initial_cells_random_markov_periodic(this_gz, mcsteps, rcell);
[fN, fnn, ~] = calc_fN(this_a0, rcell, this_gz, this_dist, this_lambda);

showlines = 0;
%h = plot_boundaries_K_Con(fN, fnn, pw, L, showlines);
interaction = [1 2];
h = plot_boundaries_K_Con_single_plot(fN, fnn, pw, L, showlines, interaction);

% Save figure
qsave = 0;
%fname_str = strrep(sprintf('Phase_boundaries_gz_%d_a0_%.1f_lambda_%.1f_%.1f_network_%d',...
%	this_gz, this_a0, this_lambda(1), this_lambda(2)), '.', 'p');
%fname_str = strrep(sprintf('Phase_boundaries_gz_%d_a0_%.1f_l_%.1f_%.1f',...
%	this_gz, this_a0, this_l(1), this_l(2)), '.', 'p');
fname_str = strrep(sprintf('Phase_boundaries_interaction_%d_%d_gz_%d_a0_%.1f_l_%.1f_%.1f',...
	interaction(1), interaction(2), this_gz, this_a0, this_l(1), this_l(2)), '.', 'p');

fname = fullfile(save_folder, 'phase_boundaries', fname_str);
save_figure(h, 5, 4, fname, '.pdf', qsave);
%% Plot Q-value vs a0
%{
idx_l = [1 1]; % fix lambda
this_l = [l1_all(idx_l(1)) l2_all(idx_l(2))];
%}
idx_lambda = [1 1];
this_lambda = [lambda1_all(1) lambda2_all(1)];

h=figure;
hold on
ls = '--';
gz_idx_all = [10 15 20 30 40 50]-4;
for idx_gz=gz_idx_all %1:numel(gz_all)    
    plot(a0_all, Q_value_all(idx_gz, :, idx_lambda(1), idx_lambda(2)), 'o--', 'LineWidth', 2);
end
xlabel('a_0');
ylabel('Q-value');
legend(sprintfc('gz = %d', gz_all(gz_idx_all)));
set(gca, 'FontSize', 32);
xlim([0 max(a0_all)]);

% Save figure
qsave = 1;
fname_str = strrep(sprintf('Q_value_vs_a0_diff_gz_lambda_%.1f_%.1f_network_%d_v3',...
	this_lambda(1), this_lambda(2), network), '.', 'p');
fname = fullfile(save_folder, fname_str);
save_figure(h, 10, 8, fname, '.pdf', qsave);
%% Plot Q-value vs gz
%{
idx_l = [1 1]; % fix lambda
this_l = [l1_all(idx_l(1)) l2_all(idx_l(2))];
%}
idx_lambda = [1 1];
this_lambda = [lambda1_all(1) lambda2_all(1)];

h=figure;
hold on
a0_idx_all = 6; %[1 5 15 40]; %[1 5 10 15 20 30 40 50]; %1:numel(a0_all); %[1 5 40] %[1 3 10];
for idx_a0=a0_idx_all 
    plot(gz_all, Q_value_all(:, idx_a0, idx_lambda(1), idx_lambda(2)), 'o--', 'LineWidth', 2)
end
xlabel('Grid size');
ylabel('max. Q-value');
set(gca, 'FontSize', 32);
legend(sprintfc('a_0 = %.1f', a0_all(a0_idx_all)), 'Location', 'ne', 'FontSize', 20);
%legend({sprintf('L=%d', L), 'L=\infty',...
%    'L=\infty, N=\infty', 'L=\infty, N=\infty, NNA'}, 'Location', 'ne');
xlim([0 50]);
ylim([0 0.05]);

% Save figure
qsave = 1;
fname_str = strrep(sprintf('Q_value_vs_gz_diff_a0_lambda_%.1f_%.1f_network_%d_v2_a0_optimal',...
	this_lambda(1), this_lambda(2), network), '.', 'p');
fname = fullfile(save_folder, fname_str);
save_figure(h, 10, 8, fname, '.pdf', qsave);

%% Plot Q-value vs gz, a0 simultaneously 
idx_lambda = [1 1];
this_lambda = [lambda1_all(1) lambda2_all(1)];

h=figure;
hold on
Q_data = squeeze(Q_value_all(6:end, :, idx_lambda(1), idx_lambda(2)));
imagesc(a0_all, gz_all(6:end), Q_data)
set(gca, 'YDir', 'normal')
xlabel('a_0');
ylabel('Grid Size');
set(gca, 'FontSize', 32);
c=colorbar;
ylim(c, [0 0.05]);
xlim([a0_all(1) a0_all(end)]);
ylim([10 gz_all(end)]);
ylabel(c, 'Q-value');

% Find maximum for each grid size and plot location (with a0 value)
Q_max_vals = max( Q_data, [], 2 );
[~, Q_max_idx] = find( Q_data == repmat(Q_max_vals, 1, numel(a0_all)) );
a0_max = a0_all(Q_max_idx);
scatter( a0_max, gz_all(6:end), 'rx' );

% Save figure
qsave = 1;
fname_str = strrep(sprintf('Q_value_vs_a0_gz_heat_map_lambda_%.1f_%.1f_network_%d_all_v1',...
	this_lambda(1), this_lambda(2), network), '.', 'p');
fname = fullfile(save_folder, fname_str);
save_figure(h, 12, 8, fname, '.pdf', qsave);
%% Plot Q-value vs lambda or l
idx_a0 = 1; % fix a0
idx_gz = 1; % fix gz
this_a0 = a0_all(idx_a0);
this_gz = gz_all(idx_gz);

h=figure;
hold on
%Q_data = squeeze(Q_value_all(idx_gz, idx_a0, :, 1:end));
%imagesc(l2_all, l1_all(1:end), Q_data)
%xlabel('l^{(2)} = \lambda^{(2)}/a_0');
%ylabel('l^{(1)} = \lambda^{(1)}/a_0');
Q_data = squeeze(Q_value_all(idx_gz, idx_a0, :, 1:end))';
imagesc(l1_all, l2_all, Q_data)
xlabel('l^{(1)} = \lambda^{(1)}/a_0');
ylabel('l^{(2)} = \lambda^{(2)}/a_0');
set(gca, 'YDir', 'normal')
set(gca, 'FontSize', 32);
%title(sprintf('gz = %d, a0 = %.2f', this_gz, this_a0));
title(sprintf('gz = %d', this_gz));
c=colorbar;
ylim([0 5]);
ylim(c, [0 0.05]);
ylabel(c, 'Q-value');

% Find maximum
[i1, i2] = find(Q_data == max(Q_data(:)));
fprintf('max(Q) = %.3f, l1 = %.2f, l2 = %.2f \n', max(Q_data(:)), l1_all(i2), l2_all(i1) )
% Plot maximum
scatter(  l1_all(i2), l2_all(i1), 100, 'rx' );
text( l1_all(i2)+0.1, l2_all(i1)+0.2,  sprintf('(%.2f, %.2f)',  l1_all(i2), l2_all(i1)), 'FontSize', 16 );

% Save figure
qsave = 0;
fname_str = strrep(sprintf(...
    'Q_value_vs_lambda_gz_%d_network_%d_range1_%.2f_%.2f_range2_%.2f_%.2f_with_max',...
	this_gz, network, l1_all(1), l1_all(end), l2_all(1), l2_all(end)), '.', 'p');
fname = fullfile(save_folder, fname_str);
save_figure(h, 12, 8, fname, '.pdf', qsave);
%% Find all maxima (over l1, l2) for different gz, a0
[Q_max_temp, idx_2_temp] = max(Q_value_all, [], 4);
[Q_max_all, idx_1] = max(Q_max_temp, [], 3);
idx_2 = diag(squeeze(idx_2_temp(:,1,idx_1)));
disp( [l1_all(idx_1); l2_all(idx_2)] )

h=figure;
plot( gz_all, Q_max_all, 'o-' );
ylim([0 0.06]);
xlabel('Grid Size');
ylabel('max Q-value');

% Find locations of maxima
h=figure;
scatter3( l2_all(idx_2), l1_all(idx_1), gz_all );
ylabel('l^{(1)}');
xlabel('l^{(2)}');
zlabel('gz');
ylim([0 l1_all(end)]);
xlim([0 l2_all(end)]);

%% Plot area fractions: example

% Plot phase boundaries
this_l = [2.2 1000];
this_lambda = this_l*this_a0;
fprintf('gz = %d, a0 = %.1f, l = [%.1f %.1f] \n',...
    this_gz, this_a0, this_l(1), this_l(2));
%}                
% -----------Calculate interaction strengths-----------------------
[fN, fnn, ~] = calc_fN(this_a0, rcell, this_gz, this_dist, this_lambda);

showlines = 1;
h = plot_boundaries_K_Con(fN, fnn, pw, L, showlines);

% Calc area fractions 
A_frac = zeros(2);
A_frac(1,1) = calc_area(L, fN, fnn, pw, 'self', 1, 0, f_MF);
A_frac(1,2) = calc_area(L, fN, fnn, pw, 'mutual', 2, 0, f_MF);
A_frac(2,1) = calc_area(L, fN, fnn, pw, 'mutual', 1, 0, f_MF);
disp(A_frac)

%% Plot area fractions: all data

A_frac = zeros(numel(l1_all), numel(l2_all), 2, 2);
for i=1:numel(l1_all)
    for j=1:numel(l2_all)
        this_l = [l1_all(i) l2_all(j)]; %[2.2 0.1];
        this_lambda = this_l*this_a0;
        fprintf('gz = %d, a0 = %.1f, l = [%.1f %.1f] \n',...
            this_gz, this_a0, this_l(1), this_l(2));
        %}                
        % -----------Calculate interaction strengths-----------------------
        [fN, fnn, ~] = calc_fN(this_a0, rcell, this_gz, this_dist, this_lambda);

        % Plot area fractions 
        A_frac(i,j,1,1) = calc_area(L, fN, fnn, pw, 'self', 1, 0, f_MF);
        A_frac(i,j,1,2) = calc_area(L, fN, fnn, pw, 'mutual', 2, 0, f_MF);
        A_frac(i,j,2,1) = calc_area(L, fN, fnn, pw, 'mutual', 1, 0, f_MF);
        %disp(A_frac)
    end
end
%%
i=1; j=2; % interaction parameters

h=figure;
hold on
A_frac_data = squeeze(A_frac(:, :, i, j))';
imagesc(l1_all, l2_all, A_frac_data)
%plot([0 6], [0.7 0.7], 'r--'); % 1<-1 and 2<-1
%plot([2.2 2.2], [0 5], 'r--'); % 1<-1 and 2<-1
xlabel('l^{(1)} = \lambda^{(1)}/a_0');
ylabel('l^{(2)} = \lambda^{(2)}/a_0');
set(gca, 'YDir', 'normal')
set(gca, 'FontSize', 32);
title(sprintf('%d \\leftarrow %d, gz = %d', i, j, this_gz));
c=colorbar;
ylim([0 5]);
%caxis([0.3 0.6]);
set(c, 'ylim', [0.3 0.6]);
ylabel(c, 'Area fraction');

% Save figure
qsave = 0;
fname_str = strrep(sprintf(...
    'Area_fraction_%d%d_vs_lambda_gz_%d_network_%d_range1_%.2f_%.2f_range2_%.2f_%.2f',...
	i,j,this_gz, network, l1_all(1), l1_all(end), l2_all(1), l2_all(end)), '.', 'p');
fname = fullfile(save_folder, fname_str);
save_figure(h, 12, 8, fname, '.pdf', qsave);

%% Area fraction "cuts" 
h11 = figure;
idx_2 = 7;
A_frac_data = squeeze(A_frac(:, idx_2, 1, 1));
plot(l1_all, A_frac_data, 'LineWidth', 2);
xlabel('l^{(1)} = \lambda^{(1)}/a_0');
set(gca, 'YDir', 'normal')
set(gca, 'FontSize', 32);
xlim([0 6]);
ylim([0 0.6]);

h12 = figure;
idx_1 = 7;
A_frac_data = squeeze(A_frac(idx_1, :, 1, 2));
plot(l2_all, A_frac_data, 'LineWidth', 2);
xlabel('l^{(2)} = \lambda^{(2)}/a_0');
set(gca, 'YDir', 'normal')
set(gca, 'FontSize', 32);
xlim([0 5]);
ylim([0 0.6]);

h21 = figure;
idx_2 = 7;
A_frac_data = squeeze(A_frac(:, idx_2, 2, 1));
plot(l1_all, A_frac_data, 'LineWidth', 2);
xlabel('l^{(1)} = \lambda^{(1)}/a_0');
set(gca, 'YDir', 'normal')
set(gca, 'FontSize', 32);
xlim([0 6]);
ylim([0 0.6]);

% Save figures
qsave = 1;
fname_str = strrep(sprintf(...
    'Area_fraction_cuts_%d%d_vs_lambda_gz_%d_network_%d_range_%.2f_%.2f',...
	1, 1, this_gz, network, l1_all(1), l1_all(end)), '.', 'p');
fname = fullfile(save_folder, fname_str);
save_figure(h11, 10, 8, fname, '.pdf', qsave);

fname_str = strrep(sprintf(...
    'Area_fraction_cuts_%d%d_vs_lambda_gz_%d_network_%d_range_%.2f_%.2f',...
	1, 2, this_gz, network, l2_all(1), l2_all(end)), '.', 'p');
fname = fullfile(save_folder, fname_str);
save_figure(h12, 10, 8, fname, '.pdf', qsave);

fname_str = strrep(sprintf(...
    'Area_fraction_cuts_%d%d_vs_lambda_gz_%d_network_%d_range_%.2f_%.2f',...
	2, 1, this_gz, network, l1_all(1), l1_all(end)), '.', 'p');
fname = fullfile(save_folder, fname_str);
save_figure(h21, 10, 8, fname, '.pdf', qsave);
%% Functions
function [Q_value] = calc_Q_value(L, fN, fnn, pw, N_inf, f_MF)
    % N_inf = 1: N->infty limit
    % mean-field contribution
    if ~N_inf
        f_MF = fN'-6*fnn;
    end
    
    % calculate intersection points and define integrals
    syms Con % for integrals
    % ---self interaction---
    % intersection points 
    Cmax_self = (L - 1 - 4.*fnn - f_MF.*(1-pw))./(2.*fnn+f_MF.*pw);
    Cmin_self = (L - 1 - 6.*fnn - f_MF.*(1-pw))./(f_MF.*pw);
    
    % Integrals
    I_self = cell(2, 1);
    I_self{1} = @(Con) 2.*fnn.*(Con - 1);
    I_self{2} = @(Con) (L - 1 - 6.*fnn - f_MF.*(Con.*pw + 1-pw));
    I_self{3} = @(Con) 0;

    % ---mutual interaction---
    % intersection points 
    Cmax_mutual = (L - 0 - 2.*fnn - f_MF.*(1-pw))./(1 + 4.*fnn+f_MF.*pw);
    Cmin_mutual = (L - 1 - 6.*fnn - f_MF.*(1-pw))./(2*fnn + f_MF.*pw);

    % Integrals
    I_mutual = cell(2, 1);
    I_mutual{1} = @(Con) (1+2.*fnn).*(Con - 1);
    I_mutual{2} = @(Con) (L - 1 - 2.*Con.*fnn - 4.*fnn - f_MF.*(Con.*pw + 1-pw));
    I_mutual{3} = @(Con) 0;

    %% Interactions 1 <- 1 and 2 <- 1
    % determine integral bounds
    [Con_bounds, K_bounds] = get_all_bounds(Cmax_self(1), Cmin_self(1), Cmax_mutual(1), Cmin_mutual(1), L);
    Con_bound_vals = [1 Con_bounds];
    
    %disp('Con_bound_vals = ');
    %disp(Con_bound_vals); %-------------
    % disp([Ca_self(1), Cb_self(1), Ca_mutual(1), Cb_mutual(1)]);
    
    % calculate integral
    mol = 1; % molecule 1
    I_sum = 0; % integral sum value
    syms x
    for i=1:numel(Con_bounds)
        % calculate integral over bounds
        integrand_term1 = I_self{K_bounds(1, i)}(x);
        integrand_term2 = I_mutual{K_bounds(2, i)}(x);
        if ~(integrand_term1(mol)==0 && integrand_term2(mol)==0)
            this_I = double( int(integrand_term1(mol)*integrand_term2(mol), Con_bound_vals(i), Con_bound_vals(i+1)) );
            I_sum = I_sum + this_I;
        end
    end
    frac_mol_1 = I_sum/(L-1)^3;
    %% Interaction 1 <- 2 (new method)
    Ca = Cmax_mutual(2);
    Cb = Cmin_mutual(2);
    mol = 2;
    
    %disp('Ca Cb = ');
    %disp([Ca Cb]); %-------------
    
    if Ca<L && Cb<L
        %disp('case A');
        integrand_1 = I_mutual{1}(x);
        I1 =  double( int(integrand_1(mol), 1, Ca) );
        integrand_2 = I_mutual{2}(x);
        I2 =  double( int(integrand_2(mol), Ca, Cb) );
        frac_mol_2 = (I1 + I2)/(L-1)^2;
    elseif Ca<L && Cb>L
        %disp('case B');
        integrand_1 = I_mutual{1}(x);
        I1 =  double( int(integrand_1(mol), 1, Ca) );
        integrand_2 = I_mutual{2}(x);
        I2 =  double( int(integrand_2(mol), Ca, L) );
        frac_mol_2 = (I1 + I2)/(L-1)^2;
    else
        %disp('case C');
        integrand_1 = I_mutual{1}(x);
        I1 = double( int(integrand_1(mol), 1, L) );
        frac_mol_2 = I1/(L-1)^2;
    end

    % Final result
    disp(frac_mol_1)
    disp(frac_mol_2)
    Q_value = frac_mol_1*frac_mol_2;
end

function [A, C1, C2] = calc_area(L, fN, fnn, pw, interaction, molecule, N_inf, f_MF)
    % interaction: 'self' -> self (1); 'mutual' -> mutual (2)
    % molecule: 1 or 2    
    % N_inf = 1: N->infty limit
    
    % convert string input to numeric
    int = find( cellfun(@(x) ~isempty(x), regexp( interaction, {'self', 'mutual'}, 'match' )), 1);
    
    % specify molecule 
    fN = fN(molecule);
    fnn = fnn(molecule);
    pw = pw(molecule);
    
    % mean-field approximation?
    if ~N_inf % if not N-> infty
        f_MF = (fN'-6.*fnn); % restore original expression
    end
    
    % calculate intersection points and define integrals
    syms x x1 x2 % for integrals
    switch int
        case 1 % self
            % intersection points 
            C1 = (L - 1 - 4.*fnn - f_MF.*(1-pw))./(2.*fnn+f_MF.*pw);
            C2 = (L - 1 - 6.*fnn - f_MF.*(1-pw))./(f_MF.*pw);

            % Integrals
            I1 = @(x) (x.^2 - 2*x - 1).*fnn;
            I2 = @(x1, x2) L*(x2-x1) - f_MF.*pw/2.*(x2.^2-x1.^2) -...
                (1 + 6*fnn + f_MF.*(1-pw)).*(x2-x1);
        case 2 % mutual
            % intersection points 
            C1 = (L - 0 - 2.*fnn - f_MF.*(1-pw))./(1 + 4.*fnn+f_MF.*pw);
            C2 = (L - 1 - 6.*fnn - f_MF.*(1-pw))./(2*fnn + f_MF.*pw);

            % Integrals
            I1 = @(x) (x.^2 - 2*x + 1)/2.*(1+2*fnn); % integral 1
            I2 = @(x1, x2) L*(x2-x1) - (1 + 4*fnn + f_MF.*(1-pw)).*(x2-x1) -...
                (fnn + f_MF.*pw/2).*(x2.^2-x1.^2); % integral 2 -> check calculation
    end

    % Boundary cases
    if C1<L && C2<L
        %disp('case A');
        A = I1(C1) + I2(C1, C2);
    elseif C1<L && C2>L
        %disp('case B');
        A = I1(C1) + I2(C1, L);
        %disp(C1);
        %disp(I1(C1));
        %disp(I2(C1, L));
    else
        %disp('case C');
        A = I1(L);
    end
    
    A = A/(L-1)^2;
end

function [Con_bounds, K_bounds] = get_all_bounds(C1a, C1b, C2a, C2b, L)
    % C1a from Kmax = L for molecule 1
    % C1b from Kmin = L for molecule 1 
    % C2a from Kmax = L for molecule 2 
    % C2b from Kmin = L for molecule 2 
    
    % returns integral_bounds: 2 x n matrix, n<=5
    % specifies boundaries of each of the n regions for both molecules (see get_bounds)
    Con_bounds = sort([C1a C1b C2a C2b L]);
    Con_bounds = Con_bounds(Con_bounds <= L);

    K_bounds = zeros( 2, numel(Con_bounds));
    for i=1:numel(Con_bounds)
        upper_bound = Con_bounds(i);
        K_bounds(1, i) = get_bounds(C1a, C1b, upper_bound);
        K_bounds(2, i) = get_bounds(C2a, C2b, upper_bound);
    end
end
%squeeze( integral_bounds(1, :, :) )
%squeeze( integral_bounds(2, :, :) )

function bounds = get_bounds(Ca, Cb, b1)
    % 3 levels: convert to numbers    
    if Ca >= b1 && Cb >= b1
        bounds = 1; %case 1: [Kmin Kmax];
    elseif Ca <= b1 && Cb >= b1
        bounds = 2; %case 2: [Kmin L];
    else
        bounds = 3; %case 3:[L L];
    end
end

function [fN, fnn, fnnn] = calc_fN(a0, rcell, gz, dist, lambda)
    % calculate fN
    Rcell = a0*rcell;
    idx = gz + round(gz/2); % pick cell not at corner of grid
    dist_vec = a0*dist(idx,:);
    r = dist_vec(dist_vec>0);
    fN = zeros(2,1);
    fN(1) = sum(sinh(Rcell./lambda(1))*exp((Rcell-r)./lambda(1)).*(lambda(1)./r));
    fN(2) = sum(sinh(Rcell./lambda(2))*exp((Rcell-r)./lambda(2)).*(lambda(2)./r));
    
    % calculate fnn
    fnn = zeros(1, 2);
    rnn = a0;
    fnn(1) = sinh(Rcell./lambda(1))*exp((Rcell-rnn)./lambda(1)).*(lambda(1)./rnn) ; % calculate signaling strength
    fnn(2) = sinh(Rcell./lambda(2))*exp((Rcell-rnn)./lambda(2)).*(lambda(2)./rnn) ; % calculate signaling strength
    
    fnnn = zeros(2,2);
    rnnn = [sqrt(3); 2].*a0;
    fnnn(1,:) = sinh(Rcell./lambda(1))*exp((Rcell-rnnn)./lambda(1)).*(lambda(1)./rnnn);
    fnnn(2,:) = sinh(Rcell./lambda(2))*exp((Rcell-rnnn)./lambda(2)).*(lambda(2)./rnnn);
end

function h = plot_boundaries_K_Con(fN, fnn, pw, L, showlines)
    Con1_all = linspace(1, L, 100);
    Con2_all = linspace(1, L, 100);
    YMF1_all = (fN(1) - 6*fnn(1))*(Con1_all*pw(1) + (1-pw(1)));
    YMF2_all = (fN(2) - 6*fnn(2))*(Con2_all*pw(2) + (1-pw(2)));

    % bounds
    K11_all_upb = 1 + 2*Con1_all*fnn(1) + 4*fnn(1) + YMF1_all;
    K11_all_lwb = 1 + 6*fnn(1) + YMF1_all;
    K12_all_upb = Con2_all + (4*Con2_all + 2)*fnn(2) + YMF2_all;
    K12_all_lwb = 1 + (2*Con2_all + 4)*fnn(2) + YMF2_all;
    K21_all_upb = Con1_all + 4*Con1_all*fnn(1) + 2*fnn(1) + 4*fnn(1) + YMF1_all;
    K21_all_lwb = 1 + 2*Con1_all*fnn(1) + 4*fnn(1) + YMF1_all;

    % intersection points 
    C1_self = (L - 1 - 4.*fnn - (fN' - 6*fnn).*(1-pw))./(2.*fnn+(fN'-6.*fnn).*pw);
    C2_self = (L - 1 - 6.*fnn - (fN' - 6*fnn).*(1-pw))./((fN'-6.*fnn).*pw);

    C1_mutual = (L - 0 - 2.*fnn - (fN' - 6*fnn).*(1-pw))./(1 + 4.*fnn+(fN'-6.*fnn).*pw);
    C2_mutual = (L - 1 - 6.*fnn - (fN' - 6*fnn).*(1-pw))./(2*fnn + (fN'-6.*fnn).*pw);

    h=figure;
    subplot(2,2,1);
    hold on
    plot(K11_all_upb, Con1_all, '-', 'Color', 'r', 'LineWidth', 2 )
    plot(K11_all_lwb, Con1_all, '--', 'Color', 'r', 'LineWidth', 2 )
    if showlines
        plot([1 L], [C1_self(1) C1_self(1)], '--', 'Color', 'b');
        plot([1 L], [C2_self(1) C2_self(1)], '--', 'Color', 'b');
    end
    xlim([0 L]);
    ylim([0 L]);
    title('1 \leftarrow 1');
    xlabel('K^{(11)}');
    ylabel('C_{ON}^{(1)}');
    set(gca, 'FontSize', 20);

    subplot(2,2,2);
    hold on
    plot(K12_all_upb, Con2_all, '-', 'Color', 'r', 'LineWidth', 2 )
    plot(K12_all_lwb, Con2_all, '--', 'Color', 'r', 'LineWidth', 2 )
    if showlines
        plot([1 L], [C1_mutual(2) C1_mutual(2)], '--', 'Color', 'b');
        plot([1 L], [C2_mutual(2) C2_mutual(2)], '--', 'Color', 'b');
    end
    xlim([0 L]);
    ylim([0 L]);
    title('1 \leftarrow 2');
    xlabel('K^{(12)}');
    ylabel('C_{ON}^{(2)}');
    set(gca, 'FontSize', 20);

    subplot(2,2,3);
    hold on
    plot(K21_all_upb, Con2_all, '-', 'Color', 'r', 'LineWidth', 2 )
    plot(K21_all_lwb, Con2_all, '--', 'Color', 'r', 'LineWidth', 2 )
    if showlines
        plot([1 L], [C1_mutual(1) C1_mutual(1)], '--', 'Color', 'b');
        plot([1 L], [C2_mutual(1) C2_mutual(1)], '--', 'Color', 'b');
    end
    xlabel('K^{(21)}');
    ylabel('C_{ON}^{(1)}');
    xlim([0 L]);
    ylim([0 L]);
    title('2 \leftarrow 1');
    set(gca, 'FontSize', 20);
end

function h = plot_boundaries_K_Con_single_plot(fN, fnn, pw, L, showlines, interaction)
    Con_all = linspace(1, L, 100);
    YMF1_all = (fN(1) - 6*fnn(1))*(Con_all*pw(1) + (1-pw(1)));
    YMF2_all = (fN(2) - 6*fnn(2))*(Con_all*pw(2) + (1-pw(2)));

    if all(interaction==[1 1])
        K11_all_upb = 1 + 2*Con_all*fnn(1) + 4*fnn(1) + YMF1_all;
        K11_all_lwb = 1 + 6*fnn(1) + YMF1_all;
        C1_self = (L - 1 - 4.*fnn - (fN' - 6*fnn).*(1-pw))./(2.*fnn+(fN'-6.*fnn).*pw);
        C2_self = (L - 1 - 6.*fnn - (fN' - 6*fnn).*(1-pw))./((fN'-6.*fnn).*pw);

        h=figure;
        hold on
        plot(K11_all_upb, Con_all, '-', 'Color', 'r', 'LineWidth', 2 )
        plot(K11_all_lwb, Con_all, '--', 'Color', 'r', 'LineWidth', 2 )
        if showlines
            plot([1 L], [C1_self(1) C1_self(1)], '--', 'Color', 'b');
            plot([1 L], [C2_self(1) C2_self(1)], '--', 'Color', 'b');
        end
        xlim([0 L]);
        ylim([0 L]);
        title('1 \leftarrow 1');
        xlabel('K^{(11)}');
        ylabel('C_{ON}^{(1)}');
        set(gca, 'FontSize', 20);
    elseif all(interaction==[1 2])
        K12_all_upb = Con_all + (4*Con_all + 2)*fnn(2) + YMF2_all;
        K12_all_lwb = 1 + (2*Con_all + 4)*fnn(2) + YMF2_all;
        C1_mutual = (L - 0 - 2.*fnn - (fN' - 6*fnn).*(1-pw))./(1 + 4.*fnn+(fN'-6.*fnn).*pw);
        C2_mutual = (L - 1 - 6.*fnn - (fN' - 6*fnn).*(1-pw))./(2*fnn + (fN'-6.*fnn).*pw);
    
        h=figure;
        hold on
        plot(K12_all_upb, Con_all, '-', 'Color', 'r', 'LineWidth', 2 )
        plot(K12_all_lwb, Con_all, '--', 'Color', 'r', 'LineWidth', 2 )
        if showlines
            plot([1 L], [C1_mutual(2) C1_mutual(2)], '--', 'Color', 'b');
            plot([1 L], [C2_mutual(2) C2_mutual(2)], '--', 'Color', 'b');
        end
        xlim([0 L]);
        ylim([0 L]);
        title('1 \leftarrow 2');
        xlabel('K^{(12)}');
        ylabel('C_{ON}^{(2)}');
        set(gca, 'FontSize', 20);
    elseif all(interaction==[2 1])
        K21_all_upb = Con_all + 4*Con_all*fnn(1) + 2*fnn(1) + 4*fnn(1) + YMF1_all;
        K21_all_lwb = 1 + 2*Con_all*fnn(1) + 4*fnn(1) + YMF1_all;
        C1_mutual = (L - 0 - 2.*fnn - (fN' - 6*fnn).*(1-pw))./(1 + 4.*fnn+(fN'-6.*fnn).*pw);
        C2_mutual = (L - 1 - 6.*fnn - (fN' - 6*fnn).*(1-pw))./(2*fnn + (fN'-6.*fnn).*pw);

        h=figure;
        hold on
        plot(K21_all_upb, Con_all, '-', 'Color', 'r', 'LineWidth', 2 )
        plot(K21_all_lwb, Con_all, '--', 'Color', 'r', 'LineWidth', 2 )
        if showlines
            plot([1 L], [C1_mutual(1) C1_mutual(1)], '--', 'Color', 'b');
            plot([1 L], [C2_mutual(1) C2_mutual(1)], '--', 'Color', 'b');
        end
        xlabel('K^{(21)}');
        ylabel('C_{ON}^{(1)}');
        xlim([0 L]);
        ylim([0 L]);
        title('2 \leftarrow 1');
        set(gca, 'FontSize', 20);
    else
        error('Invalid interaction');
    end
end
%% Previously defined integrals (wrong) 
%{
I_mutual{1} = @(x1, x2) ((x2.^2./2 - x2) - (x1.^2./2 - x1))*(1+2*this_fnn); % integral for case 1: Kmin <= K <= Kmax
I_mutual{2} = @(x1, x2) L*(x2-x1) - ( (this_fnn + f_MF.*this_pw/2).*(x2.^2-x1.^2) -...
    (1 + 4*this_fnn + f_MF.*(1-this_pw)).*(x2-x1) ) ; % integral for case 2: Kmin <= K <= L
    
I_self{1} = @(x1, x2) ((x2.^2 - x2) - (x1.^2 - 2*x1)).*this_fnn; % integral for case 1: Kmin <= K <= Kmax
I_self{2} = @(x1, x2) L*(x2-x1) - f_MF.*this_fnn/2.*(x2.^2-x1.^2) -...
    (1 + 6*this_fnn + f_MF.*(1-this_fnn)).*(x2-x1); % integral for case 2: Kmin <= K <= L
%}