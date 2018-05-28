% This script calculates the map of p_in p_eq using exact simulation, for
% trajectories that differ from each other only by a single spin
close all
clear variables
warning off

% Parameters of the system
gridsize = 15;
N = gridsize^2;
a0 = 1.5;
Rcell = 0.2*a0;
K = 1;
Con = 1;
%flip = 1; % number of cells to flip
%--IMPORTANT---
final = 0;
% final = 0: sensitivity to initial conditions
% final = 1: stability of final conditions
%--------------

% use hexagonal lattice
[dist, pos] = init_dist_hex(gridsize, gridsize);

% Calculate the signaling strength
dist_vec = dist(1,:);
r = a0*dist_vec(dist_vec>0); % exclude self influence
fN = sum(Rcell*sum(exp(Rcell-r)./r)); % calculate signaling strength
%kon = N*(K-fN-Con)/(Con-1)/fN;
%koff = N*(K-fN-1)/(Con-1)/fN;

%---loop over # flips-------------------
%%
for flip = 1:1
%flip = 1;
% filename of the saved data
fname_str = strrep(sprintf('vs_pin_N%d_Con_%d_K_%d_a0_%.2f_flip_%d', ...
    N, Con, K, a0, flip), '.', 'p');
if final
    folder = fullfile(pwd, 'figures', 'spin_flip_final');
else
    folder = fullfile(pwd, 'figures', 'spin_flip_ini');
end
%% calculate the map
%{
[peqcount, dHcount, t_av, I_av, dH_mean, dH_std] = ...
    spin_flip_multiple_flips(dist, Con, K, a0, Rcell, flip, final);
p = (0:N)./N;
peqprob = transpose(peqcount./repmat(sum(peqcount,2),1,N+1));
dHprob = transpose(dHcount./repmat(sum(dHcount,2),1,N+1));

% save mat file
qsave = 1;
if qsave
    v=1;
    fname = fullfile(folder, 'data',...
        strcat(fname_str,'-v', int2str(v), '.mat'));
    while exist(strcat(fname, '.pdf')) == 2
        v=v+1;
        fname = fullfile(folder, 'data',...
            strcat(fname_str,'-v', int2str(v), '.mat'));
    end
    save(fname);
end
%}
%% Load previous file 
%
path = fullfile(folder, 'N225_a0_1p50_K_1to30_Con_1to30', 'data');
version = int2str(1);
%fname_in = fullfile(path, strcat(fname_str, '-v', version, '.mat'));
fname_in = fullfile(path, strcat(fname_str, '-v1.mat') );
load(fname_in);
close all
%}
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
%line = 58/N;
%plot([line line], [0 1], 'r--');

% Organize and save
save_fig = 0; % save figure? 0:no, 1: yes
if save_fig > 0
    v=1;
    fname = fullfile(folder, 'vs_pin',...
    	strcat(fname_str, '-v', int2str(v), '_pin-pout_map'));
    while exist(strcat(fname, '.pdf'), 'file') == 2
        v=v+1;
        fname = fullfile(folder, 'vs_pin',...
            strcat(fname_str, '-v', int2str(v), '_pin-pout_map'));
    end
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
xlabel('$$p_{in}$$', 'FontSize', 24)
ylabel('$$ \langle t_{eq} \rangle$$', 'FontSize', 24)

% Save
save_fig = 0; % save figure? 0:no, 1: yes
if save_fig > 0
    v=1;
    fname = fullfile(folder, 'vs_pin',...
    	strcat(fname_str,'-v', int2str(v), '_teq_av'));
    while exist(strcat(fname, '.pdf'), 'file') == 2
        v=v+1;
        fname = fullfile(folder, 'vs_pin',...
            strcat(fname_str,'-v', int2str(v), '_teq_av'));
    end
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
%{
h4 = figure(4);
errorbar(p, dH_mean, dH_std, 'LineWidth', 2);
set(gca,'FontSize', 24)
xlabel('p');
ylabel('$$d_H$$', 'Interpreter', 'latex');
%}
h4 = figure(4);
plot(p, dH_std, 'LineWidth', 2);
set(gca,'FontSize', 24)
xlabel('$$p$$');
ylabel('$$\sigma(d_H)$$', 'Interpreter', 'latex');

% Save
save_fig = 0; % save figure? 0:no, 1: yes
if save_fig > 0
    %fname = fullfile(pwd, 'figures', 'single_spin_flip', 'vs_pin',...
    %	strcat(fname_str,'-v', int2str(version), '_Hamming_avg'));
    v=1;
    fname = fullfile(folder, 'vs_pin',...
    	strcat(fname_str,'-v', int2str(v), '_std_lattice_flip'));
    while exist(strcat(fname, '.pdf'), 'file') == 2
        v=v+1;
        fname = fullfile(folder, 'vs_pin',...
            strcat(fname_str,'-v', int2str(v), '_std_lattice_flip'));
    end
    save_figure_pdf(h4, 10, 8, fname);
    save_figure_eps(h4, 10, 8, fname);
end
%% Plot heat map
h5 = figure(5);
hold on
im_fig = imagesc(p, p, dHprob);
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
ylabel('$$d_H/N$$', 'FontSize', 24)
xlim([0 1]);
ylim([0 1]);

%save 
save_fig = 0; % save figure? 0:no, 1: yes
if save_fig > 0
    v=1;
    fname = fullfile(folder, 'vs_pin',...
    	strcat(fname_str,'-v', int2str(v), '_Hamming_map'));
    while exist(strcat(fname, '.pdf'), 'file') == 2
        v=v+1;
        fname = fullfile(folder, 'vs_pin',...
            strcat(fname_str,'-v', int2str(v), '_Hamming_map'));
    end
    save_figure_pdf(h5, 10, 8, fname);
    save_figure_eps(h5, 10, 8, fname);
end
%% Final I of trajectories
h6 = figure(6);
plot(p, I_av, 'r-o')
set(gca,'FontSize', 20)
title(sprintf('$$N = %d, K = %.1f, S_{ON} = %.1f, a0 = %.1f, R = %.1f$$', ...
    N, K, Con, a0, Rcell),'FontSize', 18, 'Interpreter', 'latex')
xlabel('$$p_{in}$$', 'FontSize', 24)
ylabel('$$ \langle t_{eq} \rangle$$', 'FontSize', 24)

%save 
save_fig = 0; % save figure? 0: no, 1: yes
if save_fig > 0
    v=1;
    fname = fullfile(folder, 'vs_pin',...
    	strcat(fname_str,'-v', int2str(v), '_Iav'));
    while exist(strcat(fname, '.pdf'), 'file') == 2
        v=v+1;
        fname = fullfile(folder, 'vs_pin',...
            strcat(fname_str,'-v', int2str(v), '_Iav'));
    end
    save_figure_pdf(h6, 10, 8, fname);
    save_figure_eps(h6, 10, 8, fname);
end

end