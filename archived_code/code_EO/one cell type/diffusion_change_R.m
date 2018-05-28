% Compare the stationary state diffusion solution for changes in cell radius

clear variables
close all

R = linspace(0.1, 2, 1000);
eta = 2;
D = 0.1;
i = [10, 300, 800];
gamma = 1;
lambda = sqrt(D./gamma);
cR = eta* gamma/4/pi./R./lambda./(lambda + R);

h1 = figure(1);
plot(R, cR, 'r', 'LineWidth', 2)
hold on
scatter(R(i), cR(i), 'bx')
hold off
set(gca, 'FontSize', 15)
ylabel('c_R (a.u.)', 'FontSize', 18)
xlabel('R (a.u.)', 'FontSize', 18)

set(h1,'Units','Inches');
set(h1, 'Position', [0 0 7 6 ])
pos = get(h1,'Position');
set(h1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
fig_file = fullfile(pwd, 'figures', 'cR_R');
print(h1, fig_file,'-dpdf','-r0')


h2 = figure(2);
hold on
for k = 1:numel(i)
    r = linspace(R(i(k)), 5, 1000);
    tmp = cR(i(k))*R(i(k))*exp((R(i(k))-r)/lambda)./r;
    plot(r, tmp, 'LineWidth', 2)
end
hold off
box on
set(gca, 'yscale', 'log', 'FontSize', 15)
ylabel('c(r) (a.u.)', 'FontSize', 18)
xlabel('r (a.u.)', 'FontSize', 18)

set(h2,'Units','Inches');
set(h2, 'Position', [0 0 7 6 ])
pos = get(h2,'Position');
set(h2,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
fig_file = fullfile(pwd, 'figures', 'cR_R_cofr');
print(h2, fig_file,'-dpdf','-r0')
