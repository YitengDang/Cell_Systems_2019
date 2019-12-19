% This script calculates the map of p_in p_eq using exact simulation, for
% trajectories that differ from each other only by a single spin
close all
clear variables
%warning off
set(0, 'defaulttextinterpreter', 'tex');
%%
% Parameters of the system
gridsize = 11;
N = gridsize^2;
a0 = 5.4;
Rcell = 0.2*a0;
K = 8;
Con = 16;
hill = 2;
prec = 8;
%noise = 0.1; 
noiselist = 0.1; %10.^(-2:0.2:0);
sampmethod = 'in_out_binary_v2';

% use hexagonal lattice
[dist, pos] = init_dist_hex(gridsize, gridsize);

% Calculate the signaling strength
dist_vec = dist(1,:);
r = a0*dist_vec(dist_vec>0); % exclude self influence
fN = sum(Rcell*sum(exp(Rcell-r)./r)); % calculate signaling strength

% Find uniform lattice fixed points
fp = zeros(3, 1);
x0 = [0.03 0.2 0.65]; %estimates based on previous graph
%hfunc = @update_function_uniform;
hfunc2 = @(x) update_function_uniform(x, hill, Con, K, fN) - x; % update function(x) = x
for idx=1:3
    fp(idx) = fzero(hfunc2, x0(idx));
end
%% calculate the map
for i=1:numel(noiselist)
    ini_noise = noiselist(i);
    
    % Calculate map
    % (1) pin-pout map (continuous)
    %{
    [peqcount, dHcount, t_av, I_av, dH_mean, dH_std] = ...
        count_eq_parallel_finiteHill_monte_carlo_ini_noise(dist, Con, K, a0, Rcell, hill, prec, noise);
    label = 'pin_pout';
    %}
    % (2) fin-fout map (mapped to binary cells)
    n_smpl = 1000;
    [peqcount, dHcount, t_av, I_av, dH_mean, dH_std] =...
        count_eq_parallel_finiteHill_fin_fout_ini_noise(dist, Con, K, a0, Rcell,...
        hill, prec, fp, n_smpl, ini_noise);
    label = 'fin_fout';
    
    p = (0:N)./N;
    peqprob = transpose(peqcount./repmat(sum(peqcount,2),1,N+1));
    dHprob = transpose(dHcount./repmat(sum(dHcount,2),1,N+1));

    %% save mat file
    qsave = 1;
    if qsave
        version = 1;
        % filename of the saved data
        fname_str = strrep(sprintf('%s_perturb_ini_N%d_Con%d_K%d_a0_%.2f_hill%.2f_%s_ini_noise%.3f', ...
            label, gridsize^2, Con, K, a0, hill, sampmethod, ini_noise), '.', 'p');
        %save_folder = fullfile('H:\My Documents\Multicellular automaton\figures\finite_Hill\perturb_ini',...
        %    'vs_noise_fin_fout_binary');
        save_folder = fullfile('H:\My Documents\Multicellular automaton\figures\finite_Hill\perturb_ini',...
            'data_binary_constraint_v2');
        %fname = fullfile(pwd, 'figures', 'sensitivity_ini', 'data',...
        %    strcat(fname_str,'-v', int2str(version), '.mat'));
        fname = fullfile(save_folder, strcat(fname_str,'-v', int2str(version), '.mat'));
        while exist(fname, 'file') == 2
            version=version+1;
            %fname = fullfile(pwd, 'figures', 'sensitivity_ini', 'data',...
            %    strcat(fname_str,'-v', int2str(version), '.mat'));
            fname = fullfile(save_folder, strcat(fname_str,'-v', int2str(version), '.mat'));
        end
        save(fname);
    end
end
%}
%% Load data
%{
folder = 'H:\My Documents\Multicellular automaton\figures\finite_Hill\perturb_ini\data';
%fname_str = 'vs_ini_noise_N121_Con_16_K_8_a0_5p40_p0_0p30_noise10exp-2p0_0p0-v3'; % continuous cells
fname_str = 'vs_pin_N121_Con16_K8_a0_5p40_hill2p00_montecarlo_noise0p100-v2';
load(fullfile(folder, fname_str))
%}
%{
path = 'H:\My Documents\Multicellular automaton\figures\single_spin_flip\data';
version = int2str(1);
fname_in = fullfile(path, strcat(fname_str, '-v', version, '.mat'));
load(fname_in);
%}
%% pin-pout map
% plot the map
h1 = figure(1);
box on
hold on
im_fig = imagesc(p, p, peqprob);
% set title and font
%title(sprintf('$$N = %d, K = %.1f, S_{ON} = %.1f, a0 = %.1f, R = %.1f$$', ...
%    N, K, Con, a0, Rcell),'FontSize', 18, 'Interpreter', 'latex')
set(gca,'Ydir','normal','FontSize', 32)
% set invisible parts where count is zero
set(im_fig, 'AlphaData', peqcount' > 0);
% set colorbar and labels
colormap('winter')
c = colorbar;
c.Label.String = 'Probability';
%xlabel('p_{in}', 'FontSize', 32)
%ylabel('p_{out}', 'FontSize', 32)
xlabel('f_{in}');
ylabel('f_{out}');
xlim([0 1]);
ylim([-0.01 1.01]);

% Plot line p = 1/2 - 4B/fN
%B = (Con+1)/2*(1+fN) - K;
%line = 1/2 - B/(4*fN);
%plot([line line], [0 1], 'r--');

% Organize and save
save_fig = 1; % save figure? 0:no, 1: yes
if save_fig > 0
    %fname = fullfile(pwd, 'sensitivity_ini', 'vs_pin',...
    % 	strcat(fname_str,'-v', int2str(version), '_pin-pout_map'));
    %folder = 'H:\My Documents\Thesis\Sensitivity to initial conditions\Fig2';
    folder = 'H:\My Documents\Multicellular automaton\figures\finite_Hill\perturb_ini\binary_v2_example';
    fname_str = sprintf('N121_Con16_K8_a0_5p40_hill2p00_montecarlo_noise0p100_n_smpl_%d-v1_pin-pout_map', n_smpl);
    fname = fullfile(folder, fname_str);
    save_figure_pdf(h1, 10, 8, fname);
    %save_figure_eps(h1, 10, 8, fname);
end
%%
%{
% Plot the average number of steps it takes to reach equilibrium
h2 = figure(2);
plot(p, t_av, 'r-o')
set(gca,'FontSize', 20)
title(sprintf('$$N = %d, K = %.1f, S_{ON} = %.1f, a0 = %.1f, R = %.1f$$', ...
    N, K, Con, a0, Rcell),'FontSize', 18, 'Interpreter', 'latex')
xlabel('$$k_{in}$$', 'FontSize', 24)
ylabel('Average # steps for eq.', 'FontSize', 24)

% Save
save_fig = 0; % save figure? 0:no, 1: yes
if save_fig > 0
    fname = fullfile(pwd, 'sensitivity_initial_cond_continuum', 'vs_pin',...
    	strcat(fname_str,'-v', int2str(version), '_teq_av'));
    save_figure_pdf(h2, 10, 8, fname);
    save_figure_eps(h2, 10, 8, fname);
end
%}
%{
h3 = figure(3);
% calculate the analytical formula and plot
[~,omegak] = entropy_eq_sphere(dist_vec, Con, K, a0, Rcell);
plot(p , log(omegak), 'LineWidth', 1.5)
set(gca,'FontSize', 20)
xlabel('p');
ylabel('$$\Omega_k$$', 'Interpreter', 'latex');
%}
%% (Hamming) distance between final configurations
% of unflipped and flipped configurations
%{
% Plot average and std
h4 = figure(4);
errorbar(p, dH_mean, dH_std, 'LineWidth', 2);
set(gca,'FontSize', 24)
xlabel('$$p$$');
ylabel('$$\langle d_H \rangle$$', 'Interpreter', 'latex');

% Save
save_fig = 0; % save figure? 0:no, 1: yes
if save_fig > 0
    fname = fullfile(pwd, 'sensitivity_ini', 'vs_pin',...
    	strcat(fname_str,'-v', int2str(version), '_Hamming_avg'));
    save_figure_pdf(h4, 10, 8, fname);
    save_figure_eps(h4, 10, 8, fname);
end
%}
%% Plot heat map of (Hamming) distance
h5 = figure(5);
box on
hold on
im_fig = imagesc(p, p, dHprob);
% set title and font
%title(sprintf('$$N = %d, K = %.1f, S_{ON} = %.1f, a0 = %.1f, R = %.1f$$', ...
%    N, K, Con, a0, Rcell),'FontSize', 18, 'Interpreter', 'latex')
set(gca,'Ydir','normal','FontSize', 32)
% set invisible parts where count is zero
set(im_fig, 'AlphaData', dHcount' > 0);
% set colorbar and labels
colormap('winter')
c = colorbar;
c.Label.String = 'Probability';
xlabel('f_{in}', 'FontSize', 32)
ylabel('d_H/N', 'FontSize', 32)
xlim([0 1]);
ylim([-0.01 1.01]);

%save 
save_fig = 1; % save figure? 0:no, 1: yes
if save_fig > 0
    %fname = fullfile(pwd, 'sensitivity_ini', 'vs_pin',...
    % 	strcat(fname_str,'-v', int2str(version), '_Hamming_map'));
    %folder = 'H:\My Documents\Thesis\Sensitivity to initial conditions\Fig2';
    folder = 'H:\My Documents\Multicellular automaton\figures\finite_Hill\perturb_ini\binary_v2_example';
    fname_str = sprintf('N121_Con16_K8_a0_5p40_hill2p00_montecarlo_noise0p100_n_smpl_%d-v1_Hamming_map', n_smpl);
    fname = fullfile(folder, fname_str);
    save_figure_pdf(h5, 10, 8, fname);
    %save_figure_eps(h5, 10, 8, fname);
end

%% Final I of trajectories
h6 = figure(6);
plot(p, I_av, 'r-o')
set(gca,'FontSize', 20)
title(sprintf('$$N = %d, K = %.1f, C_{ON} = %.1f, a0 = %.1f, R = %.1f$$', ...
    N, K, Con, a0, Rcell),'FontSize', 18, 'Interpreter', 'latex')
xlabel('$$p_{in}$$', 'FontSize', 24)
ylabel('$$\langle I_{eq} \rangle $$', 'FontSize', 24)
%
%save
save_fig = 0; % save figure? 0:no, 1: yes
if save_fig > 0
    %fname_str = strcat(fname_str,'-v', int2str(version), '_Iav');
    fname_str = '11';
    fname = fullfile(pwd, 'sensitivity_ini', 'vs_pin',...
    	fname_str);
    fname = 'H:\My Documents\Multicellular automaton\finite_Hill\sensitivity_ini\1';
    save_figure_pdf(h6, 10, 8, fname);
    save_figure_eps(h6, 10, 8, fname);
end
%}