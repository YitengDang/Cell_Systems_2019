% Plot the function n choose k for several values of N and save in PDF

close all
clear variables
warning off

h = figure(1);
ax1 = subplot(3, 1, 1);
N = 81;
p = linspace(0, 1, N+1);
data = zeros(size(p));
for i = 1:numel(p)
    data(i) = nchoosek(N, round(p(i)*N));
end
plot(p, data, 'k', 'Linewidth', 2)
set(gca, 'FontSize', 16)
ylim([0 max(data)])

ax2 = subplot(3, 1, 2);
N = 121;
p = linspace(0, 1, N+1);
data = zeros(size(p));
for i = 1:numel(p)
    data(i) = nchoosek(N, round(p(i)*N));
end
plot(p, data, 'k', 'Linewidth', 2)
set(gca, 'FontSize', 16)
ylim([0 max(data)])
ylabel('Number of possible combinations', 'FontSize', 18)

ax3 = subplot(3, 1, 3);
N = 441;
p = linspace(0, 1, N+1);
data = zeros(size(p));
for i = 1:numel(p)
    data(i) = nchoosek(N, round(p(i)*N));
end
plot(p, data, 'k', 'Linewidth', 2)
set(gca, 'FontSize', 16)
ylim([0 max(data)])
xlabel('p = n/N', 'FontSize', 18)

set(h,'Units','Inches');
set(h, 'Position', [0 0 8 15 ])
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
out_file = fullfile(pwd, 'figures', 'nchoosek');
print(h, out_file,'-dpdf','-r0')
