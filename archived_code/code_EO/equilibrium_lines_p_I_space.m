% Trying to find the condition in p and I that gets to equilibrium. It led
% me not far (Eduardo 27-7-2016)
close all
clear variables

% Parameters of the system
gridsize = 21;
N = gridsize^2;
a0 = 0.5;
Rcell = 0.2*a0;
save_fig = 0;

% use hexagonal lattice
[pos,ex,ey] = init_cellpos_hex(gridsize,gridsize);
dist = dist_mat(pos,gridsize,gridsize,ex,ey);
dist_vec = dist(1,:);
r = a0*dist_vec(dist_vec>0); % exclude self influence

K = 10;
Son = 5;

% Get the signalling length
fN = sum(sinh(Rcell)*sum(exp(Rcell-r)./r)); % calculate signaling strength sphere

p = 0:0.01:1;

I1 = @(p) fN + (fN+1-K)./(Son-1)./p;
I2 = @(p) (K-Son-fN)./(1-p)./(Son-1) - p.*fN./(1-p);

figure(1)
plot(p, I1(p), '-r');
hold on
plot(p, I2(p), '-g');
hold off
ylim([-1 1])