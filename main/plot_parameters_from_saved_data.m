% Plot the macroscopic parameters of saved data
clear variables
close all

% Load data
[fname, path, ~] = uigetfile(fullfile(pwd,'data','dynamics'));
load(fullfile(path,fname))

% Get the data
for i=1:t
    cells = cells_hist{i};
    p(i) = sum(cells)/numel(cells);
    I(i) = moranI(cells, a0*dist);
end

subplot(3,1,1)
plot(0:50, p(1:51), '-r', 'LineWidth', 2)
%xlabel('Time steps')
set(gca, 'FontSize', 14)
ylabel('p', 'FontSize', 20)
set(gca,'XTick',[]);

subplot(3,1,2)
plot(0:50, I(1:51), '-r', 'LineWidth', 2)
%xlabel('Time steps')
set(gca, 'FontSize', 14)
ylabel('I', 'FontSize', 20)
ylim([0 0.6])
set(gca,'XTick',[]);

subplot(3,1,3)
plot(0:50, -mom(1:51)/N, '-r', 'LineWidth', 2)
set(gca, 'FontSize', 14)
xlabel('Time steps', 'FontSize', 20)
ylabel('h', 'FontSize', 20)
ylim([-7 0])

save_figure_svg(gcf, 6, 9, fullfile(pwd,'figures','orscillation_parameters.svg'))
