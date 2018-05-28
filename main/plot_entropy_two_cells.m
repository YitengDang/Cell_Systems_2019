clear variable
close all

data_path = '~/Dropbox/Matlab codes/data_twocelltypes_entropy';

[fnamein, pathin, ~] = uigetfile(fullfile(data_path));
load(fullfile(pathin,fnamein))

dist_vec = dist_vec*a0;
r = dist_vec(dist_vec>0); % exclude self influence
fN = sum(sinh(Rcell)*sum(exp(Rcell-r)./r)); % calculate signaling strength
C = (N1*Son1 + (N-N1)*Son2)/N;
    % Plot the result
    h1 = figure(1);
plot(K2,omega, 'r-', 'Linewidth', 1)
hold on
plot(K2sim,omegasim, 'ko')
% Uncomment to test limit lines
% yl = ylim;
% K = 1 + fN*C;
% plot([K K], yl, 'r')
% K = Son2 + fN*C;
% plot([K K], yl, 'k')
% K = 1 + fN;
% plot([K K], yl, 'b')
% K = Son2 + fN;
% plot([K K], yl, 'g')
hold off
legend({'Analytic','Simulation'})
xlabel('K_2')
ylabel('Entropy')
set(gca, 'FontSize', 16)

[~, filename, ~] = fileparts(fnamein);

out_file = fullfile(pwd, 'figures', 'twocelltypes', filename);
save_figure_svg(h1, 10, 6, out_file);