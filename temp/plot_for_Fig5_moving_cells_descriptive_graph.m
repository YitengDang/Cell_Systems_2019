%% Plot for use as figure describing how moving cells / randomized lattice works
x_cells = [-2 -1 0 1 2];
y_cells = [0 0 0 0 0];

% 
x = [-2.5 -2:.01:2 2.5];
norm1 = normpdf(x,0, 0.1);
norm2 = normpdf(x,0, 0.5);

h=figure;
hold on
scatter(x_cells, y_cells, 200, 'filled');
p1 = plot( x, norm1, 'LineWidth', 2 );
p2 = plot( x, norm2, 'LineWidth', 2 );
xlabel('Position', 'Interpreter', 'none');
ylabel('P(x_{new})', 'Interpreter', 'tex');
%plot( x, norm3 );
set(gca, 'XTick', [], 'YTick', []);
title({'x_{i,new} = x_{i,old} + dx_i', 'dx_i \sim N(0,\sigma_D)'}, 'Interpreter', 'tex');
xlim([-2.5 2.5]);
ylim([-0.1 4]);
%legend(p1, '\sigma_D = 0.1');
%legend(p2, '\sigma_D = 0.5');
legend([p1, p2], sprintfc('\\sigma_D = %.1f', [0.1 0.5]));
set(gca, 'FontSize', 24);

qsave = 1;
save_path_fig = 'H:\My Documents\Multicellular automaton\paper_2_draft\figures\originals\Fig6-extensions';
fname_str = 'moving_cells_sigma_D_distribution_add_text';
fname = fullfile(save_path_fig, fname_str);
save_figure(h, 10, 8, fname, '.pdf', qsave);

%% Failed attempt
sigma_D_all = [0.001 0.003 0.005 0.01 0.03 0.05 0.1 0.5];
dimensions = 2*[2.205 2.204]; %width, height
x1 = [74.775 80.255];
x2 = [80.282 80.255];
x3 = [77.529 75.532];

dx12 = sqrt(sum((x1-x2).^2));
dx13 = sqrt(sum((x1-x3).^2));
a0 = 2*mean([dx12 dx13]);

sigma_D_pixels = sigma_D_all.*a0;
circle_widths = mean(dimensions) + sigma_D_pixels*2;

circle_width = circle_widths(end);
new_pos = [75.05 79.153];
new_pos_center = new_pos + dimensions/2;
circle_pos = new_pos_center - circle_width/2;
