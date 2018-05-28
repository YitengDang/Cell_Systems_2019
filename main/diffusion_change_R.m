% Compare the stationary state diffusion solution for changes in cell radius.
% This evaluates if the form of the function is changed by a change in R

clear variables
close all

% Parameters of the run
R = linspace(0.1, 2, 1000); % Parameters of cell radius to test
eta = 2; % secretion rate
D = 0.1; % diffusion constant
i = [10, 300, 800]; % index of the R to plot as points and analysed
gamma = 1; % degradation rate
lambda = sqrt(D./gamma); % diffusion length
cR = eta* gamma/4/pi./R./lambda./(lambda + R); % Constant to calculate secretion function

h1 = figure(1);
% plot the constant vs R
plot(R, cR, 'r', 'LineWidth', 2)
hold on
% Plot the scatter points
scatter(R(i), cR(i), 'bx')
hold off
% Set fonts and labels
set(gca, 'FontSize', 15)
ylabel('c_R (a.u.)', 'FontSize', 18)
xlabel('R (a.u.)', 'FontSize', 18)

% Save figure as pdf
set(h1,'Units','Inches');
set(h1, 'Position', [0 0 7 6 ])
pos = get(h1,'Position');
set(h1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
fig_file = fullfile(pwd, 'figures', 'cR_R'); % filename
print(h1, fig_file,'-dpdf','-r0')

% For each scatter, plot the form of the function in distance
h2 = figure(2);
hold on
for k = 1:numel(i)
    r = linspace(R(i(k)), 5, 1000); % The distance vector changes with R
    tmp = cR(i(k))*R(i(k))*exp((R(i(k))-r)/lambda)./r;
    % Plot the function shape
    plot(r, tmp, 'LineWidth', 2)
end
hold off
box on
% Set scales and labels
set(gca, 'yscale', 'log', 'FontSize', 15)
ylabel('c(r) (a.u.)', 'FontSize', 18)
xlabel('r (a.u.)', 'FontSize', 18)

% Save as pdf
set(h2,'Units','Inches');
set(h2, 'Position', [0 0 7 6 ])
pos = get(h2,'Position');
set(h2,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
fig_file = fullfile(pwd, 'figures', 'cR_R_cofr'); % filename
print(h2, fig_file,'-dpdf','-r0')
