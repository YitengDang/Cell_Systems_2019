close all

% plot the Multicellular Hamiltonian map

% Parameters of the system
gridsize = 21;
N = gridsize^2;
a0 = 0.5;
Rcell = 0.2*a0;
iniON = 10;
% parameters
Son = 8;
K = 16;

% Initialize parameters
if gridsize == 50
    load('dist5050.mat')
else
    [pos,ex,ey] = init_cellpos_hex(gridsize,gridsize);
    dist = dist_mat(pos,gridsize,gridsize,ex,ey);
end

dist_vec = a0*dist(1,:);
r = dist_vec(dist_vec>0); % exclude self influence
fN = sum(sinh(Rcell)*sum(exp(Rcell-r)./r)); % calculate signaling strength

E = @(p, I) -0.5*(Son-1)*(1 + 4*fN.*p.*(1-p).*I + fN*(2*p-1).^2) ...
    -(2*p-1).*(0.5*(Son+1)*(1+fN) - K);

u = @(p, I) 2*(Son-1)*fN*(2*p-1).*(1-I)+(Son+1)*(1+fN)-2*K;
v = @(p, I) 2*(Son-1)*fN*p.*(1-p);

piv = (0:N)/N;
Iv = -0.05:0.1:1;

[p_i,I] = meshgrid(piv, Iv);
contourf(piv, Iv, E(p_i, I))