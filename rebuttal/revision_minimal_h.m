%% Calculates the minimal hamiltonian in the restricted region (by Imax)
% ABANDONED

clear all
close all

% Parameters
gridsize = 11;
N = gridsize^2;
a0=0.5;
Rcell = 0.2*a0;
Con = 8;
K = 16;

% calculate fN
[pos,ex,ey] = init_cellpos_hex(gridsize,gridsize);
dist = dist_mat(pos,gridsize,gridsize,ex,ey);
dist_vec = a0*dist(1,:);
r = dist_vec(dist_vec>0); % exclude self influence
fN = sum(sinh(Rcell)*sum(exp(Rcell-r)./r)); % calculate signaling strength

% Hamiltonian
E = @(p, I) -0.5*(Con-1)*(1 + 4*fN.*p.*(1-p).*I + fN*(2*p-1).^2) ...
    -(2*p-1).*(0.5*(Con+1)*(1+fN) - K);
piv = (0:N)/N;
Iv = -0.05:0.05:1;
[p_i,Imesh] = meshgrid(piv, Iv);
h1 = figure();
contourf(piv, Iv, E(p_i, Imesh),'LineStyle', 'none')
colormap('summer')
c = colorbar;
ylim([-0.05 1]);

% Imax 
figure(h1);
hold on
fa0 = sinh(Rcell)*exp(Rcell-a0)/a0;
Imax = @(p) (6 - 4./sqrt(p*N) - 6./(p*N) - 6*p)*fa0./(1-p)/fN;
Imax2 = @(p) (6*p - 4./sqrt((1-p)*N) - 6./((1-p)*N))*fa0./p/fN;
%plot(piv, Imax(piv), 'b')
plot(piv, max(Imax(piv), Imax2(piv)), '--b', 'Linewidth', 1.5)