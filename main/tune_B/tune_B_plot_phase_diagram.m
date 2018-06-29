% Plot phase diagram with B protocol on heatmap of B
clear variables
close all
clc
%%
% Options
protocol = 6; %name of the protocol
protid = num2str(protocol);
save = 1;

% Parameters
gridsize = 11;
N = gridsize^2;
a0 = 0.5;
Rcell = 0.2*a0;

% Get phenotype regions
[dist, pos] = init_dist_hex(gridsize, gridsize);
dist_vec = dist(1,:);
r = a0*dist_vec(dist_vec>0); % exclude self influence

% Get the signalling length
fN = sum(sinh(Rcell)*sum(exp(Rcell-r)./r)); % calculate signaling strength sphere

% Specify protocol
% Protocol 1, 3
%{
Son_ini = 5;
tsteps = 10;
K_ini = (Son_ini+1)/2*(1+fN);
K = linspace(K_ini, K_ini, tsteps+1);
if find(protocol == [1 5], 1)
    Son = linspace(5, 15, tsteps+1);
elseif protocol == 3
    load('tuneB_protocol3_K4p5to4p5_Son5to11p5_B0to5.mat');
end
%}

% Protocol 4
%{
tsteps = 10;
K_ini = 12;
Son_ini = 2*K_ini/(1+fN) - 1;
K = linspace(K_ini,2,tsteps+1);
Son = linspace(Son_ini, Son_ini, tsteps+1);
%}
% Protocol 6
%
tsteps = 10;
K_ini = 12;
Son_ini = 2*K_ini/(1+fN) - 1;
K = linspace(K_ini,4,tsteps+1);
Son = linspace(Son_ini, Son_ini, tsteps+1);
%}
%%
% Plot the result
set(0,'defaulttextinterpreter','latex');
h = figure(1);
hold on
K_x = linspace(1,20,1000);
Son_y= linspace(1,30,1000);
[X,Y] = meshgrid(K_x, Son_y);
B = (Y+1)/2*(1+fN) - X;
imagesc(K_x, Son_y, B);
set(gca, 'ydir', 'normal', 'FontSize', 20) 
c = colorbar;
colormap('winter');

plot([fN+1 fN+1], [1 Son_y(end)], 'w', 'LineWidth', 1.5) % ON region
plot(K_x, K_x./(1+fN), 'k', 'LineWidth', 1.5) % OFF region
%if (K(end/2)-fN) < (K(end/2)-1)/fN && (K(end/2)-1)/fN > Son(end/2)
    plot(K_x, K_x-fN, 'r', 'LineWidth', 1.5) % autonomous 1
    plot(K_x, (K_x-1)/fN, 'g', 'LineWidth', 1.5) % autonomous 2
%end
%plot(K, K-1, 'k--');
plot(K_x, ones(size(K_x))*(1+fN)/(1-fN), 'w--') % deactivation-autonomy transition at B=0
plot((Son_y+1)/2*(1+fN), Son_y, 'k--'); % B=0 
%plot(K, K, '--k', 'LineWidth', 1.5) % Son = K
% Plot protocol
plot(K, Son, 'x-', 'Color', 'yellow', 'MarkerSize', 8)
plot(K(1), Son(1), '^','Color', 'yellow', 'MarkerSize', 10)
plot(K(end), Son(end), 'o', 'Color', 'yellow', 'MarkerSize',10)

% Figure ticks
tick_K = [1 5 10 15 20];
tick_Son = [1 5 10 15 20 25 30];

% Figure markup
hold off
ylim([1 Son_y(end)])
xlim([1 K_x(end)])
%ylim(c, [-30 15]);
ylabel(c, 'B', 'FontSize', 24, 'rot', 90);
set(gca, 'xtick', tick_K, 'ytick', tick_Son)
xlabel('$$K$$', 'FontSize', 24)
ylabel('$$C_{ON}$$', 'FontSize', 24)

% Organize and save
if save > 0
    fname = strrep(sprintf('phase_diagram_K_Con_with_B_a0_%.1f_protocol_%s',...
        a0, protid), '.','p');
    set(h,'Units','Inches');
    set(h, 'Position', [0 0 10 6])
    pos = get(h,'Position');
    set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    fig_file = fullfile(pwd, 'figures', 'work_distribution', fname);
    print(h, fig_file,'-dpdf','-r0')
end
%%
disp((Son_y-1)/(Son_y+1));
disp(fN);