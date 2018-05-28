%% Simulates the PDE for one signalling molecule 
% (1) quasi-steady state approximation: one equation for c
close all 
clear all
clc
%%
% parameters
global d hill K Con A gamma
d = 0.001;
hill = Inf;
K = 5;
Con = 8;
A = 10;
gamma = 1;

% I/O
fname_str = strrep(sprintf('PDE1_D%.3f_hill%.1f_K%d_Con%d_gamma%d_tauC%d_A%d_bc_Neumann_0', ...
    d, hill, K, Con, gamma, 1, A), '.','p');
all_qsave = 0;

% solve equation
m = 0;
xmesh = linspace(0, 10, 100);
tspan = linspace(0, 20, 100);

sol = pdepe(m,@pdefun,@icfun,@bcfun,xmesh,tspan);
c = sol(:,:,1);

%% A surface plot is often a good way to study a solution.
% NB secretion rate S is a linearly increasing function of C, so plotting S
% instead of C should not give huge differences
h=figure;
%{
s=surf(xmesh,tspan,c);
s.EdgeColor = 'none';
lable = '_surf_plot_c';
%}
%
imagesc(xmesh, tspan, c);
set(gca, 'YDir', 'normal');
cbar = colorbar;
cmap = colormap();
%set(cbar, 'Ticks', [0 1]);
ylabel(cbar, 'Concentration');
label = '_cplot_c';
%}

title(sprintf('D = %.3f, hill = %.2f', d, hill));
xlabel('Distance x');
ylabel('Time t');
zlabel('Concentration c');
set(h, 'Units', 'Inches', 'Position', [0 0 8 8]);
set(gca, 'FontSize', 15);

qsave = 1;
if qsave && all_qsave
    fname = fullfile(pwd,'figures', strcat(fname_str, label));
    save_figure(h, 8, 6, fname, '.pdf');
end

%% Plot cell state (=secretion rate) in (x,t)
if hill==Inf
    s = heaviside(c-K);
else
    s = c^hill/(K^hill + c^hill);
end

h1=figure;
%s=surf(xmesh,tspan,s);
%s.EdgeColor = 'none';
imagesc(xmesh, tspan, s);
set(gca, 'YDir', 'normal');
title(sprintf('D = %.3f, hill = %.2f', d, hill));
xlabel('Distance x');
ylabel('Time t');
zlabel('Secretion rate s');
set(h1, 'Units', 'Inches', 'Position', [0 0 8 6]);
set(gca, 'FontSize', 15);
cbar = colorbar;
cmap = colormap();
colormap(cmap([1 end], :));
set(cbar, 'Ticks', [0 1]);
ylabel(cbar, 'Cell State');

qsave = 1;
if qsave && all_qsave
    fname = fullfile(pwd,'figures', strcat(fname_str, '_cplot_s'));
    save_figure(h1, 8, 6, fname, '.pdf');
end
%% Also study profile at end
h2=figure;
plot(xmesh, sol(end, :, 1), 'bx')
title(sprintf('D = %.3f, hill = %.2f', d, hill));
xlabel('Distance x');
ylabel('Concentration c');
set(gcf, 'Units', 'Inches', 'Position', [9 0 8 6]);
set(gca, 'FontSize', 15);
ylim([0 10]);

qsave = 0;
if qsave && all_qsave
    fname = fullfile(pwd,'figures', strcat(fname_str, '_t_final'));
    save_figure(h2, 8, 6, fname, '.pdf');
end

%% Also study initial profile
h3=figure;
plot(xmesh, sol(1, :, 1), 'rx')
title(sprintf('D = %.3f, hill = %.2f', d, hill));
xlabel('Distance x');
ylabel('Concentration c');
set(gcf, 'Units', 'Inches', 'Position', [9 0 8 6]);
set(gca, 'FontSize', 15);
ylim([0 10]);

qsave = 0;
if qsave && all_qsave
    fname = fullfile(pwd,'figures', strcat(fname_str, '_t_initial'));
    save_figure(h3, 8, 6, fname, '.pdf');
end
%% --------------------------------------------------------------------------

function [g,f,sigma] = pdefun(x, t, c, DcDx)
    global d hill K Con gamma
    tauC = 1;
    %d = 1;
    %gamma = 1;
    %K = 5;
    %Con = 8;
    Coff = 1;
    %hill = 2;
    
    if hill==Inf
        seq = (Con-Coff)*heaviside(c-K) + Coff;
    else
        seq = (Con-Coff)*c^hill/(K^hill + c^hill) + Coff;
    end
    
    g = tauC;
    f = d*DcDx;
    sigma = -gamma*c + seq;
end
% --------------------------------------------------------------------------

function u0 = icfun(x)
    global A
    %A = 10;
    u0 = A/2*(1+sin(pi/2*x));
end

% --------------------------------------------------------------------------

function [pl,ql,pr,qr] = bcfun(xl, cl, xr, cr, t)
    pl = 0; %cl;
    ql = 1;
    pr = 0; %cr;
    qr = 1;
end