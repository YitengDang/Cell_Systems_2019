% Compare fN with its estimative fN_est for different grid sizes
% This is calculated for several grid sizes to show how fN depends on N.
% See appendix of the thesis.
clear variables
close all
a0_vec = 5; %[0.5 1 1.5]; % Different a0 to test
colors = {'k','r','b'};
gridsize = 4:1:40; % gridsizes to test

fN = zeros(numel(gridsize), 1);
fN_est = zeros(numel(gridsize), 1);

for k = 1:numel(a0_vec)
    disp(k);
    h = figure(k);
    hold on
    a0 = a0_vec(k);
    Rcell = 0.2*a0;
    for i = 1:numel(gridsize)
        % use hexagonal lattice
        [dist, pos] = init_dist_hex(gridsize(i), gridsize(i));
        dist_vec = dist(1,:);
        r = dist_vec(dist_vec>0); % exclude self influence
        fN(i) = sum(Rcell*sum(exp(Rcell-a0*r)./(a0*r)));
        fN_est(i) = max_moranI(a0,Rcell,numel(dist_vec),1);
    end
    % Plot the comparison of estimatives
    plot(gridsize, fN, '-r', 'LineWidth', 1.5)
    plot(gridsize, fN_est, '--k', 'LineWidth', 1.5)
    % Set legends, fonts and labels
    legend({'Exact', 'Approx.'}, 'Location', 'se')
    set(gca, 'FontSize', 20)
    xlabel('Grid Size', 'FontSize', 24)
    ylabel('$$f_N$$', 'FontSize', 24)
    ylim([0.01 0.03]);
    %title(sprintf('a_0 = %.1f',a0), 'FontSize', 24)
    hold off
    %% Save figure as pdf
    outfile = fullfile(pwd, 'temp', strrep(sprintf('fN_vs_N_a0%.1f',a0), '.','p'));
    save_figure_pdf(h, 7, 6, outfile);
end