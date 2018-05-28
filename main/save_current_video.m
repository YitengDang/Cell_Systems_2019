close all

% Save data from the last visualization in video format and plot the graphs
% of the parameters

videoname = fullfile(pwd, 'videos', ...
    sprintf('N%d_n%d_nf%d_K%d_Con%d_t%d_noise%d.avi', N, iniON, sum(cells), K, Son, t, noise));
v = VideoWriter(videoname);
open(v);
writeVideo(v, F);
close(v);

h1 = figure;
plot(0:t, -mom/N, '-r', 'LineWidth', 1.5)
set(gca, 'FontSize', 20)
xlabel('Time (steps)', 'FontSize', 24)
ylabel('H/N', 'FontSize', 24)
outfile = fullfile(pwd, 'figures', ...
    sprintf('N%d_n%d_nf%d_K%d_Con%d_t%d_noise%d_H.pdf', N, iniON, sum(cells), K, Son, t, noise));
save_figure_pdf(h1, 8, 6, outfile)

h2 = figure;
plot(0:t, I, '-r', 'LineWidth', 1.5)
set(gca, 'FontSize', 20)
xlabel('Time (steps)', 'FontSize', 24)
ylabel('Spatial order I', 'FontSize', 24)
outfile = fullfile(pwd, 'figures', ...
    sprintf('N%d_n%d_nf%d_K%d_Con%d_t%d_noise%d_I.pdf', N, iniON, sum(cells), K, Son, t, noise));
save_figure_pdf(h2, 8, 6, outfile)

h3 = figure;
plot(0:t, Non/N, '-r', 'LineWidth', 1.5)
hold on
pest = normcdf(0, K-Son*(1+fN), noise);
pm = mean(Non/N);
plot([0 t], [pm pm], '--b', 'LineWidth', 1.5)
plot([0 t], [pest pest], '--k', 'LineWidth', 1.5)
hold off
ylim([0.98 1])
legend({'Sim.', 'Mean Sim.', 'Est. mean'},'Location','south')
set(gca, 'FontSize', 20)
xlabel('Time (steps)', 'FontSize', 24)
ylabel('p=N_{ON}/N', 'FontSize', 24)
outfile = fullfile(pwd, 'figures', ...
    sprintf('N%d_n%d_nf%d_K%d_Con%d_t%d_noise%d_p.pdf', N, iniON, sum(cells), K, Son, t, noise));
save_figure_pdf(h3, 8, 6, outfile)

figure(9)
[~,omegak] = entropy_eq_sphere(dist_vec, Son, K, a0, Rcell);
for k = 0:N
    omegak(k+1) = omegak(k+1)*nchoosek(N, k);
end
plot((0:N)/N , omegak, 'LineWidth', 1.5)