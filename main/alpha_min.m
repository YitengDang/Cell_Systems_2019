% Estimate the minimum noise that has influence on the system
% It assumes that the minimum noise that affect the system is such that the
% the mean concentration Yi - K is more than 1/sqrt(N) than the transition
% from ON to OFF (Yi - K = 0)

close all
clear variables
warning off

% Parameters of the system
gridsize = 15;
N = gridsize^2;
a0 = 0.5;
Rcell = 0.2*a0;
pv = [1]; % calculate for several fractions of ON cells

% parameters of the genetic circuit
%Son = 1:0.5:30;
%K = 1:0.5:20;
Son = 8:9;
K = 15:16;

% Initialize the grid of cells
[dist, pos] = init_dist_hex(gridsize, gridsize);

% The distance vector are the same for all the cells (periodic boundary)
dist_vec = dist(1,:);
r = a0*dist_vec(dist_vec>0); % exclude self influence
fN = sum(sinh(Rcell)*sum(exp(Rcell-r)./r)); % calculate the signaling strength

h1 = figure(1);
for p = pv
    % For each combination of Son and K calculate the minimum noise
    for i = 1:numel(Son)
        for j = 1:numel(K)
            Ynei = fN*(p*Son(i) + 1 - p);
            alpha(i,j) = min(abs(Ynei+Son(i)-K(j)), abs(Ynei+1-K(j)))/sqrt(N)/K(j);
        end
    end
    % Plot the data in the figure h1
    clf(h1,'reset');
    figure(h1)
    contourf(K, Son, log10(alpha), 'LineStyle', 'none')
    % Set axes, legend, colorbar, labels
    set(gca, 'ydir', 'normal')
    colormap('summer')
    c = colorbar;
    set(gca, 'FontSize', 20)
    xlabel('K', 'FontSize', 24)
    ylabel('C_{ON}', 'FontSize', 24)
    ylabel(c, '\alpha', 'FontSize', 24)
    % Save the figure
    %fname = sprintf('alphamin_a0%d_N%d_p%d.pdf', 10*a0, N, 10*p);
    %save_figure_pdf(h1, 8, 6, fullfile(pwd, 'figures', fname))
end
