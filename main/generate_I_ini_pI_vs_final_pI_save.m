% This script calculates P(p_out, I_out | p_in, I_in)
close all
clear all
%warning off

% Parameters of the system
gridsize = 15;
N = gridsize^2;
%a0 = 0.5;
%Rcell = 0.2*a0;
%K = 15;
%Con = 8;
a0all = [0.5 1.5 1.5];
Kall = [15 6 20];
Conall = [8 15 15];

%%
%for trial=1:3
trial = 3;
a0 = a0all(trial);
Rcell = 0.2*a0;
K = Kall(trial);
Con = Conall(trial);
fprintf('a0 = %.1f, K = %d, Con = %d \n', a0, K, Con);

p1 = floor(0.1*N)/N;
dp = round(0.05*N)/N;
p_all = p1:dp:(1-p1);
I_all = -0.05:0.05:0.4;

% use hexagonal lattice
[dist, pos] = init_dist_hex(gridsize, gridsize);

% Calculate the signaling strength
dist_vec = dist(1,:);
r = a0*dist_vec(dist_vec>0); % exclude self influence
fN = sum(sinh(Rcell)*sum(exp(Rcell-r)./r)); % calculate signaling strength

fprintf('Points to sample = %d \n', numel(p_all)*numel(I_all));
%% calculate the map
[p_out, I_out, t_av, errors] = count_pI_eq_parallel(dist, Con, K, a0, Rcell, p_all, I_all);

%% save mat file
fname_str = strrep(sprintf('pIin_pIout_N%d_Con_%d_K_%d_a0_%.1f', ...
    N, Con, K, a0), '.', 'p');
i=1;
fname = fullfile(pwd, 'data', 'pIin_pIout', strcat(fname_str,'-v', int2str(i), '.mat'));
while exist(fname, 'file') == 2
    i = i+1;
    fname = fullfile(pwd, 'data', 'pIin_pIout', strcat(fname_str,'-v', int2str(i), '.mat'));
end
save(fname)

%% plot expected input-output relations
% calculate expected final (p,I) for each initial (p,I)
p_out_av = mean(p_out, 3);
I_out_av = mean(I_out, 3);

% figure
h1 = figure(1);
hold on
[pmesh, Imesh] = meshgrid(p_all, I_all);
plot(pmesh, Imesh, 'ro');
for i=1:numel(pmesh)
    plot([pmesh(i) p_out_av(i)], [Imesh(i) I_out_av(i)], 'r-', 'LineWidth', 0.1);
end
plot(p_out_av, I_out_av, 'rx'); % plot final states
set(gca,'Ydir','normal','FontSize', 24)
% set colorbar and labels
xlabel('p', 'FontSize', 24)
ylabel('I', 'FontSize', 24)

% save
qsave=0;
if qsave
   fname_str = strrep(sprintf('pIin_pIout_N%d_a0_%.1f_Con%d_K%d_I_%.3f_%.3f_p_%.3f_%.3f_avg_pI_out',...
       N, a0, Con, K, I_all(1), I_all(end), p_all(1), p_all(end)), '.', 'p');
   fname = fullfile(pwd, 'figures', 'pIin_pIout', fname_str);
   save_figure_pdf(h1, 10, 8, fname);
   save_figure_eps(h1, 10, 8, fname);
end

%end