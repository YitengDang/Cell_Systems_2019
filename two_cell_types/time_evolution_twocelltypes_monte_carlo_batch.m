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
Mcomm = [1 0; 0 1]; % Communication matrix, type i reacts to type j iff M_ij=1
Mcomm = logical(Mcomm);
subfolder = sprintf('Mcomm_%d_%d_%d_%d', Mcomm(1,1), Mcomm(1,2),...
    Mcomm(2,1), Mcomm(2,2));

% Circuit parameters 
%Con = [15 8]; %column j = type j 
Coff = [1 1]; %default [1 1]
K = [5 5];

% Simulation parameters
p_ini = [0.5 0.5]; %fraction of type-1 and type-2 cells which are ON
tmax = 1000;
n_sim = 100;

% Loop parameters
Con1_all = 5:5:40;
Con2_all = 10; %5:5:40;
%[Con1_mesh, Con2_mesh] = meshgrid(Con1_all, Con2_all);

%--------------------------------------------------------------------------
%% Calculate interaction matrix
% Distance and position
[dist, pos] = init_dist_hex(gridsize, gridsize);
[f_mat, g_mat, f_std, g_std] = calc_interaction_strengths(N, frac_type1, a0, rcell, dist, Mcomm);

if sum(sum(f_std > 10^(-2)))
    disp('Warning! Interaction strengths have large standard deviation');
end

%% Dynamics with only p1, p2
p_out_all = zeros(numel(Con1_all), numel(Con2_all), n_sim, numel(p_ini) );
t_out_all = zeros(numel(Con1_all), numel(Con2_all), n_sim);
for il1=1:numel(Con1_all)
    for il2=1:numel(Con2_all)
        Con = [Con1_all(il1) Con2_all(il2)];
        fprintf('Con1 = %d, Con2 = %d \n', Con(1), Con(2));
        %disp(Con);
        for il3=1:n_sim
            %disp(il3)
            %----
            % individual simulation
            t = 0;
            [p_out, terminate] = update_monte_carlo(N1, N2, p_ini, Con, Coff, K, f_mat, g_mat);
            while ~terminate && t<=tmax
                t = t+1;
                p = p_out;
                [p_out, terminate] = update_monte_carlo(N1, N2, p, Con, Coff, K, f_mat, g_mat);
            end
            %----
            p_out_all(il1, il2, il3, :) = p_out;
            t_out_all(il1, il2, il3) = t;
        end
    end
end

%%
% Save result
qsave = 1;
if qsave
    fname_str = strrep(...
        sprintf('N1_%d_N2_%d_p1_%.2f_p2_%.2f_a0_%.2f_K1_%d_K2_%d_Con_%d_to_%d_Rcell1_%.2fa0_Rcell2_%.2fa0_nruns_%d_%s', ...
        N1, N2, p(1), p(2), a0, K(1), K(2), Con1_all(1), Con1_all(end), ...
        rcell(1), rcell(2), n_sim, subfolder), '.', 'p');
    i = 1;
    %fname = fullfile('K:\bn\hy\Shared\Yiteng\Multicellularity\data_twocelltypes\time_evolution', fname_str);
    fname = fullfile(pwd, 'data', 'time_evolution', 'montecarlo', subfolder,...
        strcat(fname_str,'-v',int2str(i),'.mat'));
    while exist(fname, 'file') == 2
        i=i+1;
      fname = fullfile(pwd, 'data', 'time_evolution', 'montecarlo', subfolder,...
          strcat(fname_str,'-v',int2str(i),'.mat'));
    end
    save(fname, 'p_out_all', 't_out_all', 'N1', 'N2', 'a0', 'rcell', 'K', ...
        'Con1_all', 'Con2_all', 'f_mat', 'g_mat', 'frac_type1', 'Mcomm')
end
%% Plot heat maps of <p1> and <p2>
h1 = figure(1);
imagesc(Con1_all, Con2_all, mean(p_out_all(:,:,:,1), 3)' )
title(sprintf('$$K_1 = %d, K_2 = %d$$', K(1), K(2)));
set(gca, 'YDir', 'normal', 'FontSize', 24);
xlabel('$$C_{ON,1}$$')
ylabel('$$C_{ON,2}$$')
c=colorbar;
ylabel(c, '$$\langle p_1(t_f) \rangle$$', 'Interpreter', 'latex')
caxis([0 1]);

qsave = 0;
if qsave
    fname_str = strrep(sprintf('p1_out_N%d_a0_%.1f_N1_%d_K1_%d_K2_%d_nruns_%d_%s',...
        N, a0, N1, K(1), K(2), n_sim, subfolder), '.', 'p');
    path = fullfile(pwd, 'figures', 'time_evolution_montecarlo', subfolder);
    fname = fullfile(path, fname_str);
    %fname = fullfile(pwd, 'data', 'time_evolution', fname_str);
    save_figure_pdf(h1, 10, 8, fname);
end

h2 = figure(2);
imagesc(Con1_all, Con2_all, mean(p_out_all(:,:,:,2), 3)' )
title(sprintf('$$K_1 = %d, K_2 = %d$$', K(1), K(2)));
set(gca, 'YDir', 'normal', 'FontSize', 24);
xlabel('$$C_{ON,1}$$')
ylabel('$$C_{ON,2}$$')
c=colorbar;
ylabel(c, '$$\langle p_2(t_f) \rangle$$', 'Interpreter', 'latex')
caxis([0 1]);

qsave = 0;
if qsave
    fname_str = strrep(sprintf('p2_out_N%d_a0_%.1f_N1_%d_K1_%d_K2_%d_nruns_%d_%s',...
        N, a0, N1, K(1), K(2), n_sim, subfolder), '.', 'p');
    path = fullfile(pwd, 'figures', 'time_evolution_montecarlo', subfolder);
    fname = fullfile(path, fname_str);
    %fname = fullfile(pwd, 'data', 'time_evolution', fname_str);
    save_figure_pdf(h2, 10, 8, fname);
end
%% Plot heat map of <t_out>
h3 = figure(3);
imagesc(Con1_all, Con2_all, mean(t_out_all(:,:,:), 3)' )
title(sprintf('$$K_1 = %d, K_2 = %d$$', K(1), K(2)));
set(gca, 'YDir', 'normal', 'FontSize', 24);
xlabel('$$C_{ON,1}$$')
ylabel('$$C_{ON,2}$$')
c=colorbar;
ylabel(c, '$$\langle t_{out} \rangle$$', 'Interpreter', 'latex')
%caxis([0 1]);

qsave = 0;
if qsave
    fname_str = strrep(sprintf('t_out_N%d_a0_%.1f_N1_%d_K1_%d_K2_%d_nruns_%d_tmax_%d_%s',...
        N, a0, N1, K(1), K(2), n_sim, tmax, subfolder), '.', 'p');
    path = fullfile(pwd, 'figures', 'time_evolution_montecarlo', subfolder);
    fname = fullfile(path, fname_str);
    %fname = fullfile(pwd, 'data', 'time_evolution', fname_str);
    save_figure_pdf(h3, 10, 8, fname);
end
%% Plot dynamics
%{
h1 = figure(1);
hold on
plot(0:t, p1_all, 'bo-');
plot(0:t, p2_all, 'ro-');
legend({'p_1', 'p_2'});
xlabel('time');
ylabel('fraction ON cells');
ylim([0 1]);
set(gca, 'FontSize', 24);
%}
%% Plot hamiltonian
%{
h2 = figure(2);
plot(0:t, h/N, 'o-');
xlabel('time');
ylabel('pseudo-energy h');
%}
%% Plot spatial index I
%{
h3 = figure(3);
plot(0:t, I, 'o-');
xlabel('time');
ylabel('spatial index I');
ylim([-0.2 1]);
%}
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
%{
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