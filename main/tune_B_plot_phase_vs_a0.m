% Plot fN change with a0 and phase borders to show how changing B (by
% changing a0) takes the system from one phase to another.
clear variables
close all

% Options
save = 1;
tick_K = [1 5 10 15 20];
tick_Son = [1 5 10 15 20 25 30];

% Parameters
gridsize = 11;
N = gridsize^2;

% Load data
%[fname, path, ~] = uigetfile(fullfile(pwd,'data'));
%load(fullfile(path,fname));

% variable a0, calculate fN
a0 = linspace(0.5, 1.5, 100);
[dist, pos] = init_dist_hex(gridsize, gridsize);
dist_vec = dist(1,:);
fN = zeros(1, 100);
for i=1:100
    r = a0(i)*dist_vec(dist_vec>0); % exclude self influence
    Rcell = 0.2*a0(i);
    fN(i) = sum(sinh(Rcell)*sum(exp(Rcell-r)./r)); % calculate signaling strength sphere
    %[out, map] = get_phenotype_map(a0, dist, Rcell, Son, K);
end

% fixed circuit parameters
Son = 5;
K = (Son+1)/2*(1+fN(end));

%% Plot B vs. a0
% Plot the result
set(0,'defaulttextinterpreter','latex');
h = figure(1);
hold on
%plot(a0,(Son+1)/2*(1+fN)-K, 'b', 'LineWidth', 1.5);
plot(a0, fN, 'r', 'LineWidth', 1.5);
% plot intersections 
plot(a0, (K-1)/Son*ones(1,100), 'g--');
plot(a0, (K-Son)*ones(1,100), 'k--');
%
ylim([-1 2.5]);
set(gca, 'FontSize', 24);
xlabel('$$a_0$$');
ylabel('Concentration');
legend({'f_N', 'K-1/C_{on}', 'K-C_{ON}'}, 'Location', 'eastoutside');
set(h,'Units','Inches');
set(h, 'Position', [0 0 10 6])
%%
% Organize and save
if save > 0
    fname = strrep(sprintf('phase_diagram_a0_K%.2f_Son%d_protocol2',...
        K, Son), '.','p');
    out_file = fullfile(pwd, 'figures', 'work_distribution', strcat(fname, '_withB'));
    save_figure_pdf(h, 10, 6, out_file);
end