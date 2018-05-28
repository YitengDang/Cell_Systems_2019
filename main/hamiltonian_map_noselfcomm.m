clear variables
close all

% Parameters of the system
gridsize = 15;
N = gridsize^2;
a0 = 0.5;
Rcell = 0.2*a0;
iniON = round(0.6*N);
% parameters
Son = 11;
K = 14;

% Initialize parameters
[pos,ex,ey] = init_cellpos_hex(gridsize,gridsize);
dist = dist_mat(pos,gridsize,gridsize,ex,ey);

dist_vec = a0*dist(1,:);
r = dist_vec(dist_vec>0); % exclude self influence
fN = sum(sinh(Rcell)*sum(exp(Rcell-r)./r)); % calculate signaling strength

h = @(p,I) -0.5*(Son-1)*(4*fN.*p.*(1-p).*I + fN*(2*p-1).^2)-(0.5*(Son+1)*(fN)-K)*(2*p-1);

piv = (0:2:N)/N;
Iv = -0.5:0.1:1;
[p_i,Imesh] = meshgrid(piv, Iv);

surf(p_i, Imesh, h(p_i,Imesh));
set(gca, 'FontSize', 15);
xlabel('p', 'FontSize', 20);
ylabel('I', 'FontSize', 20);
zlabel('h', 'FontSize', 20);