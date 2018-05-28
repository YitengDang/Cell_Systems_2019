close all
clear all
warning off

% Lattice parameters
gridsize = 25;
N = gridsize^2;
a0 = 0.5;
R = 0.2*a0;

%% Exact fN
% use hexagonal lattice
[dist, pos] = init_dist_hex(gridsize, gridsize);
dist = round(dist, 5);
dist_vec = a0*dist(1,:);
r = dist_vec(dist_vec>0); % exclude self influence
fN = sum(sinh(R)*sum(exp(R-r)./r)); % calculate signaling strength

%% approximate fN
fa0 = sinh(R)*sum(exp(R-a0)./a0);
G = 2*pi*exp(R)*sinh(R);
Acell = sqrt(3)/2*a0^2;
L = sqrt((N-7)*Acell/pi + (1.5*a0)^2);
fN2 = 6*fa0 + G/Acell*(exp(-1.5*a0) - exp(-L));
