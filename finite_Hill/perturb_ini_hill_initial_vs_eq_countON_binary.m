% This script calculates the map of p_in p_eq using exact simulation, for
% trajectories that differ from each other only by a single spin
close all
clear all
%warning off

% Parameters of the system
gridsize = 11;
N = gridsize^2;
a0 = 0.5;
Rcell = 0.2*a0;
K = 14;
Con = 16;

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
%%
set(0, 'defaulttextinterpreter', 'latex');
% plot the map
h1 = figure(1);
hold on
im_fig = imagesc(p,p,peqprob);
% set title and font
%title(sprintf('$$N = %d, K = %.1f, S_{ON} = %.1f, a0 = %.1f, R = %.1f$$', ...
%    N, K, Con, a0, Rcell),'FontSize', 18, 'Interpreter', 'latex')
set(gca,'Ydir','normal','FontSize', 24)
% set invisible parts where count is zero
set(im_fig, 'AlphaData', peqcount' > 0);
% set colorbar and labels
c = colorbar;
c.Label.String = 'Probability';
xlabel('$$p_{in}$$', 'FontSize', 24)
ylabel('$$p_{out}$$', 'FontSize', 24)
xlim([0 1]);
ylim([0 1]);

% Plot line p = 1/2 - 4B/fN
%B = (Con+1)/2*(1+fN) - K;
%line = 1/2 - B/(4*fN);
%plot([line line], [0 1], 'r--');

% Organize and save
save_fig = 1; % save figure? 0:no, 1: yes
if save_fig > 0
    fname = fullfile(pwd, 'figures', 'single_spin_flip', 'vs_pin',...
    	strcat(fname_str,'-v', int2str(version), '_pin-pout_map'));
    save_figure_pdf(h1, 10, 8, fname);
    save_figure_eps(h1, 10, 8, fname);
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
    save_figure_eps(h2, 10, 8, fname);
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
save_fig = 1; % save figure? 0:no, 1: yes
if save_fig > 0
    fname = fullfile(pwd, 'figures', 'single_spin_flip', 'vs_pin',...
    	strcat(fname_str,'-v', int2str(version), '_Hamming_avg'));
    save_figure_pdf(h4, 10, 8, fname);
    save_figure_eps(h4, 10, 8, fname);
end
%% Plot heat map
h5 = figure(5);
hold on
im_fig = imagesc(p, 0:N, dHprob);
% set title and font
%title(sprintf('$$N = %d, K = %.1f, S_{ON} = %.1f, a0 = %.1f, R = %.1f$$', ...
%    N, K, Con, a0, Rcell),'FontSize', 18, 'Interpreter', 'latex')
set(gca,'Ydir','normal','FontSize', 24)
% set invisible parts where count is zero
set(im_fig, 'AlphaData', dHcount' > 0);
% set colorbar and labels
c = colorbar;
c.Label.String = 'Probability';
xlabel('$$p_{in}$$', 'FontSize', 24)
ylabel('$$P(d_H|p_{in})$$', 'FontSize', 24)
xlim([0 1]);
ylim([0 N]);

%save 
save_fig = 1; % save figure? 0:no, 1: yes
if save_fig > 0
    fname = fullfile(pwd, 'figures', 'single_spin_flip', 'vs_pin',...
    	strcat(fname_str,'-v', int2str(version), '_Hamming_map'));
    save_figure_pdf(h5, 10, 8, fname);
    save_figure_eps(h5, 10, 8, fname);
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