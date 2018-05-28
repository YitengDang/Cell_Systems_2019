% This script calculates the map of p_in p_eq using exact simulation, for a
% system with a single burst at t=0
close all
clear all
%warning off

% Parameters of the system
gridsize = 11;
N = gridsize^2;
a0 = 0.5;
Rcell = 0.2*a0;
K = 18;
Son = 7;
B = 5000; %burst size

%subfolder in H:\My Documents\Multicellular automaton\figures\pin_pout\ to save result in
save_fig = 1; % save figure? 0:no, 1: yes
subfolder = 'deactivation burst';

% use hexagonal lattice
[dist, pos] = init_dist_hex(gridsize, gridsize);

% filename of the saved data
fname = sprintf('pin_pout_N%d_Con_%d_K_%d_a0_%d_B_%d', ...
    N, Son, round(K), 10*a0, B);

% Calculate the signaling strength
dist_vec = dist(1,:);
r = a0*dist_vec(dist_vec>0); % exclude self influence
fN = sum(Rcell*sum(exp(Rcell-r)./r)); % calculate signaling strength
kon = N*(K-fN-Son)/(Son-1)/fN;
koff = N*(K-fN-1)/(Son-1)/fN;

% calculate the map
[count, t_av] = count_eq_parallel_burst(dist, Son, K, a0, Rcell,B);
p = (0:N)./N;
prob = transpose(count./repmat(sum(count,2),1,N+1));

% save mat file
save(fullfile(pwd,'figures','pin_pout',subfolder,strcat(fname,'.mat')))
%%
% plot the map
h1 = figure(1);
im_fig = imagesc(p,p,prob);
% set title and font
title(sprintf('N = %d, K = %d, S_{ON} = %d, a0 = %.1f, R = %.1f, B=%d', ...
    N, K, Son, a0, Rcell, B),'FontSize', 18);
set(gca,'Ydir','normal','FontSize', 20)
% set invisible parts where count is zero
set(im_fig, 'AlphaData', count' > 0);
% set colorbar and labels
c = colorbar;
c.Label.String = 'Probability';
xlabel('p_{in}', 'FontSize', 24)
ylabel('p_{out}', 'FontSize', 24)

% Organize and save
if save_fig > 0
    set(h1,'Units','Inches');
    set(h1, 'Position', [0 0 10 6 ])
    pos = get(h1,'Position');
    set(h1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    out_file = fullfile(pwd, 'figures', 'pin_pout', subfolder, strcat(fname,'_map'));
    print(h1, out_file,'-dpdf','-r0')
end

% Plot the average number of steps it takes to reach equilibrium
h2 = figure(2);
plot(p, t_av, 'r-o')
set(gca,'FontSize', 20)
title(sprintf('N = %d, K = %d, S_{ON} = %d, a0 = %.1f, R = %.1f, B=%d', ...
    N, K, Son, a0, Rcell, B),'FontSize', 18)
xlabel('k_{in}', 'FontSize', 24)
ylabel('Average # steps for eq.', 'FontSize', 24)

% Optional saving:
if save_fig > 0
    set(h2,'Units','Inches');
    set(h2, 'Position', [0 0 10 6 ])
    pos = get(h2,'Position');
    set(h2,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    out_file = fullfile(pwd, 'figures', 'pin_pout', subfolder, strcat(fname,'_times'));
    print(h2, out_file,'-dpdf','-r0')
end

h3 = figure(3);
% calculate the analytical formula and plot
[~,omegak] = entropy_eq_sphere(dist_vec, Son, K, a0, Rcell);
plot((0:N)/N , log(omegak), 'LineWidth', 1.5)