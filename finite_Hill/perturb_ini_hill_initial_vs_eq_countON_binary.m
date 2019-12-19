% This script calculates the map of p_in p_eq using exact simulation, for
% trajectories that differ from each other only by a single spin
close all
clear all
%warning off
set(0, 'defaulttextinterpreter', 'tex');
%%
% Parameters of the system
gridsize = 11;
N = gridsize^2;
a0 = 0.5;
Rcell = 0.2*a0;
K = 10;
Con = 5;

% use hexagonal lattice
[dist, pos] = init_dist_hex(gridsize, gridsize);

% filename of the saved data
fname_str = sprintf('vs_pin_N%d_Con_%d_K_%d_a0_%d_binary', ...
    N, Con, round(K), 10*a0);

% Calculate the signaling strength
dist_vec = dist(1,:);
r = a0*dist_vec(dist_vec>0); % exclude self influence
fN = sum(Rcell*sum(exp(Rcell-r)./r)); % calculate signaling strength
kon = N*(K-fN-Con)/(Con-1)/fN;
koff = N*(K-fN-1)/(Con-1)/fN;

%% calculate the map
[peqcount, dHcount, t_av, I_av, dH_mean, dH_std] = spin_flip_count_eq_parallel(dist, Con, K, a0, Rcell);
p = (0:N)./N;
peqprob = transpose(peqcount./repmat(sum(peqcount,2),1,N+1));
dHprob = transpose(dHcount./repmat(sum(dHcount,2),1,N+1));

% save mat file
qsave = 1;
if qsave
    version=1;
    fname = fullfile(pwd, 'figures', 'single_spin_flip', 'data',...
        strcat(fname_str,'-v', int2str(version), '.mat'));
    while exist(fname, 'file') == 2
        version=version+1;
        fname = fullfile(pwd, 'figures', 'single_spin_flip', 'data',...
            strcat(fname_str,'-v', int2str(version), '.mat'));
    end
    save(fname);
end
%% Load saved data
% binary cells 
%
folder = 'H:\My Documents\Multicellular automaton\figures\main\spin_flip_ini\data';
fname_str = 'vs_pin_N121_Con_5_K_10_a0_5-v1'; 
load(fullfile(folder, fname_str))
p = (0:N)./N;
%}
%% Plot pin-pout map
% plot the map
h1 = figure(1);
hold on
box on
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
xlabel('p_{in}', 'FontSize', 32)
ylabel('p_{out}', 'FontSize', 32)
xlim([0 1]);
ylim([-0.01 1.01]);

% Plot line for multiple spin flips
this_p = 58/121;
plot([this_p this_p], [0 1], 'r--')
% Plot line p = 1/2 - 4B/fN
%B = (Con+1)/2*(1+fN) - K;
%line = 1/2 - B/(4*fN);
%plot([line line], [0 1], 'r--');

% Organize and save
save_fig = 0; % save figure? 0:no, 1: yes
if save_fig > 0
    %fname = fullfile(pwd, 'figures', 'single_spin_flip', 'vs_pin',...
    % 	strcat(fname_str,'-v', int2str(version), '_pin-pout_map'));
    folder = 'H:\My Documents\Thesis\Sensitivity to initial conditions\Figures_I';
    fname_str = 'I_a_2_vs_pin_N121_Con_5_K_10_a0_5-v1_pin-pout_map_v3_with_line';
    fname = fullfile(folder, fname_str);
    save_figure_pdf(h1, 10, 8, fname);
    %save_figure_eps(h1, 10, 8, fname);
end
%%
% Plot the average number of steps it takes to reach equilibrium
h2 = figure(2);
plot(p, t_av, 'r-o')
set(gca,'FontSize', 20)
title(sprintf('$$N = %d, K = %.1f, S_{ON} = %.1f, a0 = %.1f, R = %.1f$$', ...
    N, K, Con, a0, Rcell),'FontSize', 18, 'Interpreter', 'latex')
xlabel('k_{in}', 'FontSize', 24)
ylabel('Average # steps for eq.', 'FontSize', 24)

% Save
save_fig = 0; % save figure? 0:no, 1: yes
if save_fig > 0
    fname = fullfile(pwd, 'figures', 'single_spin_flip', 'vs_pin',...
    	strcat(fname_str,'-v', int2str(version), '_teq_av'));
    save_figure_pdf(h2, 10, 8, fname);
    %save_figure_eps(h2, 10, 8, fname);
end

%{
h3 = figure(3);
% calculate the analytical formula and plot
[~,omegak] = entropy_eq_sphere(dist_vec, Con, K, a0, Rcell);
plot(p , log(omegak), 'LineWidth', 1.5)
set(gca,'FontSize', 20)
xlabel('p');
ylabel('$$\Omega_k$$', 'Interpreter', 'latex');
%}
%% Hamming distance between final configurations
% between unflipped and flipped configurations

% Plot average and std
h4 = figure(4);
errorbar(p, dH_mean, dH_std, 'LineWidth', 2);
set(gca,'FontSize', 24)
xlabel('p');
ylabel('$$d_H$$', 'Interpreter', 'latex');

% Save
save_fig = 0; % save figure? 0:no, 1: yes
if save_fig > 0
    fname = fullfile(pwd, 'figures', 'single_spin_flip', 'vs_pin',...
    	strcat(fname_str,'-v', int2str(version), '_Hamming_avg'));
    save_figure_pdf(h4, 10, 8, fname);
    save_figure_eps(h4, 10, 8, fname);
end
%% Plot heat map
h5 = figure(5);
box on
hold on
im_fig = imagesc(p, p, dHprob);
plot(0:1/N:1, 1/N:1/N:1+1/N, 'r--');
plot(0:1/N:1, 1+1/N:-1/N:1/N, 'g--');

% set title and font
%title(sprintf('$$N = %d, K = %.1f, S_{ON} = %.1f, a0 = %.1f, R = %.1f$$', ...
%    N, K, Con, a0, Rcell),'FontSize', 18, 'Interpreter', 'latex')
set(gca,'Ydir','normal','FontSize', 32)
% set invisible parts where count is zero
set(im_fig, 'AlphaData', dHcount' > 0);
% set colorbar and labels
c = colorbar;
colormap('winter');
c.Label.String = 'Probability';
xlabel('p_{in}', 'FontSize', 32)
ylabel('d_H(X_{in}, X_{out})/N', 'FontSize', 32)
%ylabel('P(d_H|p_{in})', 'FontSize', 32)
xlim([0 1]);
ylim([-0.01 1.01]);

%save 
save_fig = 0; % save figure? 0:no, 1: yes
if save_fig > 0
    %fname = fullfile(pwd, 'figures', 'single_spin_flip', 'vs_pin',...
    % 	strcat(fname_str,'-v', int2str(version), '_Hamming_map'));
    folder = 'H:\My Documents\Thesis\Sensitivity to initial conditions\Fig1';
    fname_str = 'I_a_2_vs_pin_N121_Con_5_K_10_a0_5-v1_Hamming_map_v3_with_lines';
    fname = fullfile(folder, fname_str);
    
    save_figure_pdf(h5, 10, 8, fname);
    %save_figure_eps(h5, 10, 8, fname);
end

%% Final I of trajectories
h6 = figure(6);
plot(p, I_av, 'r-o')
set(gca,'FontSize', 20)
title(sprintf('$$N = %d, K = %.1f, S_{ON} = %.1f, a0 = %.1f, R = %.1f$$', ...
    N, K, Con, a0, Rcell),'FontSize', 18, 'Interpreter', 'latex')
xlabel('k_{in}', 'FontSize', 24)
ylabel('Average # steps for eq.', 'FontSize', 24)

%save 
save_fig = 0; % save figure? 0:no, 1: yes
if save_fig > 0
    fname = fullfile(pwd, 'figures', 'single_spin_flip', 'vs_pin',...
    	strcat(fname_str,'-v', int2str(version), '_Iav'));
    save_figure_pdf(h6, 10, 8, fname);
    save_figure_eps(h6, 10, 8, fname);
end