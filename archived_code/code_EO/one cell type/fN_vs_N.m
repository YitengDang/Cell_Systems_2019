% Compare fN with its estimative fN_est for different grid sizes
% This is calculated for several grid sizes
clear variables
close all
a0_vec = [0.5 1 1.5];
colors = {'k','r','b'};
gridsize = 4:2:40;

fN = zeros(numel(gridsize), 1);
fN_est = zeros(numel(gridsize), 1);

for k = 1:numel(a0_vec)
    h = figure(k);
    hold on
    a0 = a0_vec(k);
    Rcell = 0.2*a0;
    for i = 1:numel(gridsize)
        % use hexagonal lattice
        [pos,ex,ey] = init_cellpos_hex(gridsize(i),gridsize(i));
        dist = dist_mat(pos,gridsize(i),gridsize(i),ex,ey);
        dist_vec = dist(1,:);
        r = dist_vec(dist_vec>0); % exclude self influence
        fN(i) = sum(Rcell*sum(exp(Rcell-a0*r)./(a0*r)));
        fN_est(i) = max_moranI(a0,Rcell,numel(dist_vec),1);
    end
    plot(gridsize, fN, '-r', 'LineWidth', 1.5)
    plot(gridsize, fN_est, '--k', 'LineWidth', 1.5)
    legend({'Exact', 'Approx.'}, 'Location', 'southeast')
    set(gca, 'FontSize', 20)
    xlabel('Grid Size', 'FontSize', 24)
    ylabel('f_N', 'FontSize', 24)
    title(sprintf('a_0 = %.1f',a0), 'FontSize', 24)
    hold off
    outfile = fullfile(pwd,'figures',sprintf('fN_vs_N_a0%d.pdf',a0));
    save_figure_pdf(h, 7, 6, outfile);
end