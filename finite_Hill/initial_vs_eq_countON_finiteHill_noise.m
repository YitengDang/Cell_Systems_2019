% This script calculates the map of p_in p_eq using exact simulation.
close all
clear variables

% Parameters of the system
gridsize = 15;
N = gridsize^2;
%a0 = 5.6;
%Rcell = 0.2*a0;
K = 8;
Con = 16;
hill = 2;
noise = 10^(0);
a0list = [8];
tmax = 1000;
n_smpl = 100;

% method for sampling initial states; binary=ON/OFF, normal = N(p, sigma)
sampmethod = 'montecarlo'; %'binary' 'normal' 'uniform'

% use hexagonal lattice
[dist, pos] = init_dist_hex(gridsize, gridsize);

% path for saving
path = fullfile(pwd, 'data', 'pin_pout', 'noise'); %data file path
path2 = fullfile(pwd, 'figures', 'pin_pout', 'noise'); %images path
if exist(path)~=7
	disp('Save directory 1 does not exist!');
elseif exist(path2)~=7
    disp('Save directory 2 does not exist!');
end

%% calculate the map
for i=1:numel(a0list)
close all;
a0 = a0list(i);
Rcell = 0.2*a0;

% base filename for saving
fname_str = strrep(sprintf('pin_pout_noise_%.3f_N%d_Con%d_K%d_a0_%.2f_hill%.2f_tmax%d_nsmpl%d_%s', ...
    noise, N, Con, K, a0, hill, tmax, n_smpl, sampmethod), '.', 'p');

% Calculate the signaling strength
dist_vec = dist(1, :);
r = a0*dist_vec(dist_vec>0); % exclude self influence
fN = sum(Rcell*sum(exp(Rcell-r)./r)); % calculate signaling strength

% calculate map
[count, corr_av] = count_eq_parallel_hill_noise_monte_carlo_tmax(...
    dist, Con, K, Rcell, a0, hill, noise, tmax, n_smpl);
p = (0:N)./N;
prob = transpose(count./repmat(sum(count,2),1,N+1));

% save mat file
file_out = fullfile(path, strcat(fname_str,'.mat'));
save(file_out);
%}
%% Load from file
%fname = strrep(sprintf('pin_pout_noise_%.3f_N%d_Con%d_K%d_a0_%.2f_hill%.2f_%s', ...
%    noise, N, Con, K, a0, hill, sampmethod), '.', 'p');
%path = 'H:\My Documents\Multicellular automaton\finite_Hill\figures\pin_pout\data';
path = 'H:\My Documents\Multicellular automaton\finite_Hill\data\pin_pout\noise';
[fname, path, ~] = uigetfile(path);
load(fullfile(path, fname));
%load(fullfile(path,strcat(fname, '.mat')));
close all;
%}
%% plot the map
set(0,'defaulttextinterpreter', 'latex');
h1 = figure(1);
im_fig = imagesc(p,p,prob);
% set title and font
%title(sprintf('$$N = %d, K = %d, S_{ON} = %.1f, a_0 = %.1f, R = %.1f$$', ...
%    N, K, Con, a0, Rcell),'FontSize', 18)
set(gca,'Ydir','normal','FontSize', 24)
% set invisible parts where count is zero
set(im_fig, 'AlphaData', count' > 0);
% set colorbar and labels
c = colorbar;
c.Label.String = 'Probability';
xlabel('$$\langle X \rangle_{in}$$', 'FontSize', 24)
ylabel('$$\langle X \rangle_{out}$$', 'FontSize', 24)

% Organize and save
% base filename for saving
path2 = 'H:\My Documents\Multicellular automaton\temp';
fname_str = strrep(sprintf('pin_pout_noise_%.3f_N%d_Con%d_K%d_a0_%.2f_hill%.2f_tmax%d_nsmpl%d_%s', ...
    10^(-0.2), N, Con, K, a0, hill, tmax, n_smpl, sampmethod), '.', 'p');

save_fig = 0;
if save_fig > 0
    out_file = fullfile(path2, strcat(fname_str, '_map'));
    save_figure_pdf(h1, 10, 8, out_file);
    save_figure_eps(h1, 10, 8, out_file);
end

%% Solve for uniform lattice fixed points and include in plot
fp = zeros(3, 1);
x0 = [0.03 0.2 0.65]; %estimates based on previous graph
hfunc = @update_function_uniform;
hfunc2 = @(x) hfunc(x, hill, Con, K, fN);
for idx=1:3
    fp(idx) = fzero(hfunc2, x0(idx));
end
figure(h1);
hold on
plot([0 1], [fp(1) fp(1)], 'r--'); % stable fixed point
%plot([0 1], [fp(2) fp(2)], 'r--');
plot([0 1], [fp(3) fp(3)], 'r--'); % stable fixed point
plot([fp(2) fp(2)], [0 1], 'g--'); % unstable fixed point

% Organize and save
save_fig = 1;
if save_fig > 0
    out_file = fullfile(path2, strcat(fname_str, '_map_v2'));
    save_figure_pdf(h1, 10, 8, out_file);
    save_figure_eps(h1, 10, 8, out_file);
end
%%
% Plot the average number of steps it takes to reach equilibrium
%{
h2 = figure(2);
plot(p, t_av, 'r-o')
set(gca,'FontSize', 24)
%title(sprintf('$$N = %d, K = %d, S_{ON} = %d, a0 = %.1f, R = %.1f$$', ...
%    N, K, Con, a0, Rcell),'FontSize', 18)
xlabel('$$\langle X \rangle_{in}$$', 'FontSize', 24)
%ylabel('Average # steps for eq.', 'FontSize', 24)
ylabel('$$\langle t_{eq} \rangle$$');
pmax = p(t_av == max(t_av));
title(sprintf('max $$t_{eq}$$ at $$p_{in}$$ = %.2f', pmax));

save_fig = 1;
if save_fig > 0
    out_file = fullfile(path2, strcat(fname_str, '_tav'));
    save_figure_pdf(h2, 10, 8, out_file);
    save_figure_eps(h2, 10, 8, out_file);
end
%}
%%
% Plot the average corr in equilibrium
h3 = figure(3);
plot(p, corr_av, 'r-o')
set(gca,'FontSize', 24)
%title(sprintf('$$N = %d, K = %d, S_{ON} = %d, a0 = %.1f, R = %.1f$$', ...
%s    N, K, Con, a0, Rcell),'FontSize', 18)
xlabel('$$\langle X \rangle_{in}$$', 'FontSize', 24)
ylabel('$$\langle C(X_{out}) \rangle$$', 'FontSize', 24)

save_fig = 0;
if save_fig > 0
    out_file = fullfile(path2, strcat(fname_str, '_corr'));
    save_figure_pdf(h3, 10, 8, out_file);
    save_figure_eps(h3, 10, 8, out_file);
end

%close all;
end
%%
%{
h4 = figure(4);
plot(p, ntmax, 'r-o')
set(gca,'FontSize', 20)
title(sprintf('$$N = %d, K = %d, S_{ON} = %d, a0 = %.1f, R = %.1f$$', ...
    N, K, Con, a0, Rcell),'FontSize', 18)
xlabel('$$p_{in}$$', 'FontSize', 24)
ylabel('Average I in equilibrium', 'FontSize', 24)
%}
%{
h3 = figure(3);
% calculate the analytical formula and plot
[~,omegak] = entropy_eq_sphere(dist_vec, Son, K, a0, Rcell);
plot((0:N)/N , log(omegak), 'LineWidth', 1.5)
%}