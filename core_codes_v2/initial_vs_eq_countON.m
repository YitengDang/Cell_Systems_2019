% This script calculates the map of p_in p_eq using exact simulation.
close all
clear all
%warning off

% Parameters of the system
gridsize = 11;
N = gridsize^2;
a0 = 1.5;
Rcell = 0.2*a0;
K = 3;
Con = 24;
% K, Con = [6, 21]; [3, 24]

% use hexagonal lattice
[dist, pos] = init_dist_hex(gridsize, gridsize);

% filename of the saved data
save_folder = 'H:\My Documents\Thesis\Miscellaneous\3_number_of_pI_states\figures';
fname_str = sprintf('pin_pout_N%d_Con_%d_K_%d_a0_%d', ...
    N, Con, round(K), 10*a0);

% Calculate the signaling strength
dist_vec = dist(1,:);
r = a0*dist_vec(dist_vec>0); % exclude self influence
fN = sum(Rcell*sum(exp(Rcell-r)./r)); % calculate signaling strength
kon = N*(K-fN-Con)/(Con-1)/fN;
koff = N*(K-fN-1)/(Con-1)/fN;

% calculate the map
[count, t_av] = count_eq_parallel(dist, Con, K, a0, Rcell);
p = (0:N)./N;
prob = transpose(count./repmat(sum(count,2),1,N+1));

% save mat file
save(fullfile(save_folder, strcat(fname_str, '_exact_sim','.mat')))
%%
% plot the map
h1 = figure;
box on
hold on
im_fig = imagesc(p,p,prob);
% set title and font
%title(sprintf('N = %d, K = %.1f, S_{ON} = %.1f, a0 = %.1f, R = %.1f', ...
%    N, K, Con, a0, Rcell),'FontSize', 18)
set(gca,'Ydir','normal', 'FontSize', 20)
% set invisible parts where count is zero
set(im_fig, 'AlphaData', count' > 0);
% set colorbar and labels
c = colorbar;
c.Label.String = 'Probability';
xlabel('p_{in}', 'FontSize', 24)
ylabel('p_{out}', 'FontSize', 24)
colormap('viridis');
% Plot line p = 1/2 - 4B/fN
%{
B = (Con+1)/2*(1+fN) - K;
line = 1/2 - B/(4*fN);
plot([line line], [0 1], 'r--');
%}
xlim([0 1]);
ylim([0 1]);

% Organize and save
save_fig = 1; % save figure? 0:no, 1: yes
if save_fig > 0
    %{
    set(h1,'Units','Inches');
    set(h1, 'Position', [0 0 10 6 ])
    pos = get(h1,'Position');
    set(h1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    out_file = fullfile(pwd, 'rebuttal', strcat(fname,'_map'));
    print(h1, out_file,'-dpdf','-r0')
    %}
    path_out = fullfile(save_folder, fname_str);
    save_figure_pdf(h1, 10, 8, path_out)
end
%%
% Plot the average number of steps it takes to reach equilibrium
h2 = figure(2);
plot(p, t_av, 'r-o')
set(gca,'FontSize', 20)
title(sprintf('N = %d, K = %.1f, S_{ON} = %.1f, a0 = %.1f, R = %.1f', ...
    N, K, Con, a0, Rcell),'FontSize', 18)
xlabel('k_{in}', 'FontSize', 24)
ylabel('Average # steps for eq.', 'FontSize', 24)

%{
h3 = figure(3);
% calculate the analytical formula and plot
[~,omegak] = entropy_eq_sphere(dist_vec, Con, K, a0, Rcell);
plot((0:N)/N , log(omegak), 'LineWidth', 1.5)
set(gca,'FontSize', 20)
xlabel('p');
ylabel('$$\Omega_k$$', 'Interpreter', 'latex');
%}