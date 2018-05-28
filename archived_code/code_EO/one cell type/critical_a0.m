% Calculate the value of a0 that leads to fN = 1 for different cell radii
% At this value the influence of neighbours and self are balanced
% This function assumes hexagonal lattice positioning of the cells and
% periodic boundary
clear all
close all

g = [11; 21; 50]; % set grid size
c = {'r-'; 'k-'; 'b-'};

% set the range of a0 in which fN will be calculated
a0 = linspace(0.001,5,1000)';
% Set the ratios between the cell radii and a0
ratio = linspace(0.02,0.5,50)';
% Initialize the vector of critical a0
a0crit = zeros(numel(ratio), 3);

% The critical a0 for each ratio is given by the value that best
% approximates 1
figure(1)
hold on
for m = 1:numel(g)
    for k = 1:size(ratio,1)
        fn = zeros(size(a0));
        for i = 1:size(a0,1)
            gridsize = g(m); % set an odd grid size
            N = gridsize^2; % number of cells
            Rcell = ratio(k)*a0(i);
            aux = dist_hex(gridsize,a0(i)); % calculate the distances
            r = aux(aux>0); % exclude the self influence
            fN(i) = sum(sinh(Rcell)*sum(exp(Rcell-r)./r)); % calculate signaling strength
        end
        [val, idx] = min(abs(fN-1)); % find the value that best approximates 1
        a0crit(k,m) = a0(idx);
    end
    plot(ratio, a0crit(:,m), c{m}, 'LineWidth', 1.5)
    leg_str{m} = sprintf('%d x %d grid', gridsize, gridsize);
end
hold off
box on
legend(leg_str)
set(gca,'FontSize', 16)
xlabel('R_{cell}/a_0','FontSize', 18)
ylabel('a_0^{c}','FontSize', 18)
