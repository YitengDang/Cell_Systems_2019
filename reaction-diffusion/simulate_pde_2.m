%% Simulates the PDE for one signalling molecule 
% (2) coupled equations for c and s
close all 
clear all
clc
%%
% parameters
global d hill K Con A tauS gamma
d = 0.00;
hill = Inf;
K = 5;
Con = 8;
A = 10;
tauS = 10;
gamma = 1;

% I/O
fname_str = strrep(sprintf(...
    'PDE2_D%.3f_hill%.1f_tauS%.2f_K%d_Con%d_gamma%d_tauC%d_A%d_bc_c_Neumann_s_Dirichlet',...
    d, hill, tauS, K, Con, gamma, 1, A), '.','p');
all_qsave = 1;

% solve equation
xmesh = linspace(0, 10, 1000); % requires high mesh value
tspan = linspace(0, 200, 1000);

m = 0;
sol = pdepe(m,@pdefun,@icfun,@bcfun,xmesh,tspan);
c = sol(:,:,1);
s = sol(:,:,2);

%% Check solution: (1) whether it's in ss at the end
dx = xmesh(2)-xmesh(1);
c_final = sol(end,:,1);
s_final = sol(end,:,2);
if hill==Inf
    seq = (Con-1)*heaviside(c_final-K) + 1;
else
    seq = (Con-1)*c_final.^hill./(K^hill + c_final.^hill) + 1;
end
ds = seq - s_final;
dc = d*del2(c_final, dx) - gamma*c_final + s_final;

figure;
hold on
plot(dc);
plot(ds);
legend({'dc/dt(tf)', 'ds/dt(tf)'})
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
% rescale s
s1 = (s-1)/(Con-1);
% compare with value for instantaneous response
if hill==Inf
    s2 = heaviside(c-K);
else
    s2 = c^hill/(K^hill + c^hill);
end

h1=figure;
%s=surf(xmesh,tspan,s);
%s.EdgeColor = 'none';
imagesc(xmesh, tspan, s);
set(gca, 'YDir', 'normal');
title(sprintf('D = %.3f, hill = %.2f', d, hill));
xlabel('Distance x');
ylabel('Time t');
zlabel('Cell state s');
set(h1, 'Units', 'Inches', 'Position', [0 0 8 6]);
set(gca, 'FontSize', 15);
cbar = colorbar;
cmap = colormap();
colormap(cmap([1 end], :));
set(cbar, 'Ticks', [0 1]);
ylabel(cbar, 'Cell State');

%figure;
%imagesc(xmesh, tspan, s2);

qsave = 1;
if qsave && all_qsave
    fname = fullfile(pwd,'figures', strcat(fname_str, '_cplot_s'));
    save_figure(h1, 8, 6, fname, '.pdf');
end
%% Also study profile at end
h2=figure;
hold on
plot(xmesh, sol(end, :, 1), 'bx');
plot(xmesh, sol(end, :, 2), 'rx'); %, 'Color', [0.5 0 0.7])
title(sprintf('t = %d, D = %.3f, hill = %.2f', tspan(end), d, hill));
xlabel('Distance x');
ylabel('Concentration c');
set(gcf, 'Units', 'Inches', 'Position', [9 0 8 6]);
set(gca, 'FontSize', 16);
legend({'c', 's'});
%ylim([0 10]);

qsave = 1;
if qsave && all_qsave
    fname = fullfile(pwd,'figures', strcat(fname_str, '_t_final'));
    save_figure(h2, 8, 6, fname, '.pdf');
end

%% Also study initial profile
t2 = 1;
h3=figure;
hold on
plot(xmesh, sol(t2, :, 1), 'bx');
plot(xmesh, sol(t2, :, 2), 'rx'); %, 'Color', [0.7 0.5 0]);
title(sprintf('t = %d, D = %.3f, hill = %.2f', t2, d, hill));
xlabel('Distance x');
ylabel('Concentration c');
set(gcf, 'Units', 'Inches', 'Position', [9 0 8 6]);
set(gca, 'FontSize', 15);
legend({'c', 's'});
%ylim([0 10]);

qsave = 1;
if qsave && all_qsave
    fname = fullfile(pwd,'figures', strcat(fname_str, '_t_initial'));
    save_figure(h3, 8, 6, fname, '.pdf');
end
%% --------------------------------------------------------------------------

function [g,f,sigma] = pdefun(x, t, c, DcDx)
    global d hill K Con tauS gamma
    tauC = 1;
    %d = 1;
    %gamma = 1;
    %K = 5;
    %Con = 8;
    Coff = 1;
    %hill = 2;
    
    if hill==Inf
        seq = (Con-Coff)*heaviside(c(1)-K) + Coff;
    else
        seq = (Con-Coff)*c(1)^hill/(K^hill + c(1)^hill) + Coff;
    end
    
    g = [tauC; tauS];
    f = [d; 0].*DcDx;
    sigma = [-gamma*c(1)+c(2); seq-c(2)];
end
% --------------------------------------------------------------------------

function u0 = icfun(x)
    global A
    %A = 10;
    n = 2.5;
    u0 = [A/2*(1+sin(2*n*pi/A*x)); 0];
end

%--------------------------------------------------------------------------

function [pl,ql,pr,qr] = bcfun(xl, cl, xr, cr, t)
    pl = [0; cl(2)]; 
    ql = [1; 0]; 
    pr = [0; cr(2)]; 
    qr = [1; 0]; 
    %pl = [0; 0]; %[0; cl(2)];
    %ql = [1; 1]; %[1; 0];
    %pr = [0; 0]; %[0; cr(2)];
    %qr = [1; 1]; %[1; 0];
end