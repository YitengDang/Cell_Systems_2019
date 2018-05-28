close all
clear all
warning off

% Estimate the minimum noise that has influence on the system

% Parameters of the system
gridsize = 11;
N = gridsize^2;
a0 = 1.5;
Rcell = 0.2*a0;
pv = [1];

% parameters
Son = 1:0.5:30;
K = 1:0.5:20;

% Initialize parameters
[pos,ex,ey] = init_cellpos_hex(gridsize,gridsize);
dist = dist_mat(pos,gridsize,gridsize,ex,ey);

dist_vec = dist(1,:);
r = a0*dist_vec(dist_vec>0); % exclude self influence
fN = sum(sinh(Rcell)*sum(exp(Rcell-r)./r));

h1 = figure(1);
for p = pv
    for i = 1:numel(Son)
        for j = 1:numel(K)
            Ynei = fN*(p*Son(i) + 1 - p);
            alpha(i,j) = min(abs(Ynei+Son(i)-K(j)), abs(Ynei+1-K(j)))/sqrt(N)/K(j);
        end
    end
    clf(h1,'reset');
    figure(h1)
    contourf(K, Son, log10(alpha), 'LineStyle', 'none')
    %imagesc(K, Son, alpha)
    set(gca, 'ydir', 'normal')
    colormap('summer')
    c = colorbar;
    set(gca, 'FontSize', 20)
    xlabel('K', 'FontSize', 24)
    ylabel('C_{ON}', 'FontSize', 24)
    ylabel(c, 'K', 'FontSize', 24)
    %text(16, 8, '1')
    fname = sprintf('alphamin_a0%d_N%d_p%d.pdf', 10*a0, N, 10*p);
    save_figure_pdf(h1, 8, 6, fullfile(pwd, 'figures', fname))
end
