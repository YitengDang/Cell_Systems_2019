% Time evolution of a system without noise and with visualization
% Finite Hill coefficient
% v2: Compare to binary system
close all
clear all
warning off
set(0, 'defaulttextinterpreter', 'latex');

%%
% Lattice parameters
gridsize = 11;
N = gridsize^2;
a0 = 1.5;
Rcell = 0.2*a0;

% circuit parameters
K = 6;
Con_all = 1:0.1:25; %1:0.1:30;
hill = 2; % Hill coefficient

% use hexagonal lattice
[dist, pos] = init_dist_hex(gridsize, gridsize);
dist = round(dist, 5);
dist_vec = a0*dist(1,:);
r = dist_vec(dist_vec>0); % exclude self influence
fN = sum(sinh(Rcell)*sum(exp(Rcell-r)./r)); % calculate signaling strength
gN = sum(sum((sinh(Rcell)*exp(Rcell-r)./r).^2)); % calculate signaling strength

%%
% Find single cell / uniform lattice fixed points
x0 = 0:0.1:1;
fp = zeros(numel(Con_all), numel(x0));

% Check stability of fixed points
stable = zeros(numel(Con_all), numel(x0));
evals = zeros(numel(Con_all), numel(x0), N);

% update functions
%hfunc = @update_function_uniform;
%hill_f_deriv = @(x) (hill*K^hill+ ((Con-1)*x + 1)^(hill-1) )/(K^hill+ ((Con-1)*x + 1)^hill)^2;
%hill_f_deriv(fp(1))

for i = 1:numel(Con_all)
    disp(i);
    Con = Con_all(i);
    hfunc = @(x) (((Con-1)*x + 1))^hill/(K^hill+ ((Con-1)*x + 1)^hill) - x;
    hfunc2 = @(x) ((1+fN)*((Con-1)*x + 1))^hill/(K^hill+ ((1+fN)*(Con-1)*x + 1)^hill) - x;
    for j=1:numel(x0)
        fp(i, j) = fzero(hfunc2, x0(j));
        [stable(i, j), evals(i,j,:)] = stability_fun(fp(i, j)*ones(N,1), a0, Rcell, K, Con, hill, dist);
    end
end


%% Plot bifurcation diagram
colors = ['r' 'b'];

h = figure(1);
hold on
for i=1:numel(Con_all)
    %disp(Con_all(i));
    for j=1:numel(x0)
        color = colors(stable(i,j)+1);
        
        plot(Con_all(i), fp(i,j), 'o', 'Color', color);
    end
end
ylim([-0.02 1.02]);
xlim([Con_all(1) Con_all(end)]);
xlabel('$$C_{ON}$$');
ylabel('$$X$$');
set(gca, 'FontSize', 24);

%Save fig
qsave = 1;
if qsave
    folder = 'H:\My Documents\Multicellular automaton\temp';
    fname = strrep(sprintf('bifurcation_diagram_vsCon_K%d_hill%.2f_a0_%.2f_spacing_%.2f',...
        K, hill, a0, Con_all(2)-Con_all(1)), '.','p');
    out_file = fullfile(folder, fname);
    save_figure_pdf(h, 10, 8, out_file)
end

%% Plot phase diagram 
%{
% For single cell in uniform lattice
% Bistability example: K=8, Con=16, a0=1.5, hill=2
set(0,'defaulttextinterpreter', 'latex');
h=figure(1);
hold on
%S = linspace(1+fN,Con*(1+fN),200);
%fS = 1 + (Con-1)*S.^hill./(K^hill+S.^hill);
%plot(linspace(0, 1, 200), (S-1)/(Con-1), 'b-', 'LineWidth', 1.5);
%plot(linspace(0, 1, 200), (fS-1)/(Con-1), 'r-', 'LineWidth', 1.5);
X = linspace(0, 1, 200);
Y = (Con-1)*X + 1;
plot(X, X, 'b-', 'LineWidth', 1.5);
plot(X, Y.^hill./(K^hill+Y.^hill), 'r-', 'LineWidth', 1.5);
%xlabel('Sensed concentration $$Y_i$$', 'FontSize', 20);
%ylabel('Secreted concentration $$f(Y_i)$$', 'FontSize', 20);
%xlabel('$$X(t)$$');
%ylabel('$$X(t+1)$$');
xlabel('X');
ylabel('f(X)');
%legend({'X_i', 'f(X_i)'}, 'FontSize', 20, 'Location', 'nw');
set(gca, 'FontSize', 24);
xlim([0 1]);
ylim([0 1]);

% Plot together with fixed points
h;
plot([fp(1) fp(1)], [0 1], 'k--');
plot([fp(2) fp(2)], [0 1], 'k--');
plot([fp(3) fp(3)], [0 1], 'k--');

%Save fig
qsave = 0;
if qsave
    fname = strrep(sprintf('uniform_lattice_map_Con%d_K%d_hill%.2f_a0_%.2f', Con, K, hill, a0), '.','p');
    out_file = fullfile(pwd, fname);
    %out_file = fullfile(pwd, 'figures', 'finite_Hill_autonomy_transition', fname);
    save_figure_pdf(h, 10, 8, out_file)
end
%}