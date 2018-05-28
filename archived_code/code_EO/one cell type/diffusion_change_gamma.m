% Compare the stationary state diffusion solution for changes in
% degradation
clear variables
close all

R = 0.2;
eta = 2;
D = 0.1;
i = [200, 500, 800];
gamma = linspace(0.1, 5, 1000);
lambda = sqrt(D./gamma);
cR = eta* gamma/4/pi/R./lambda./(lambda + R);

h1 = figure(1);
plot(gamma, cR, 'r', 'LineWidth', 2)
hold on
scatter(gamma(i), cR(i), 'bx')
hold off
set(gca, 'FontSize', 15)
ylabel('c_R (a.u.)', 'FontSize', 18)
xlabel('\gamma (a.u.)', 'FontSize', 18)

set(h1,'Units','Inches');
set(h1, 'Position', [0 0 7 6 ])
pos = get(h1,'Position');
set(h1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
fig_file = fullfile(pwd, 'figures', 'cR_gamma');
print(h1, fig_file,'-dpdf','-r0')

r = linspace(R, 5*R, 1000);
h2 = figure(2);
for k = 1:numel(i)
    out(k,:) = cR(i(k))*R*exp((R-r)/lambda(i(k)))./r;
end    
plot(r, out, 'LineWidth', 2)
set(gca, 'yscale', 'log', 'FontSize', 15)
ylabel('c(r) (a.u.)', 'FontSize', 18)
xlabel('r (a.u.)', 'FontSize', 18)

set(h2,'Units','Inches');
set(h2, 'Position', [0 0 7 6 ])
pos = get(h2,'Position');
set(h2,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
fig_file = fullfile(pwd, 'figures', 'cR_gamma_cofr');
print(h2, fig_file,'-dpdf','-r0')
