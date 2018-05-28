% Plots the limiting cases regions in color

clear variables
close all

% Parameters of the system
gridsize = 11;
N = gridsize^2;
a0 = 0.5;
Rcell = 0.2*a0;

% Get the values for which to make the map
Son_vec = linspace(1,30,1000);
K_vec = linspace(0,20,1000);
[K, Son] = meshgrid(K_vec, Son_vec);

% Get the signalling length
fN = sum(Rcell*sum(exp(Rcell-r)./r)); % calculate signaling strength sphere

% Make 4 limiting regions as boolean matrices
R1 = 1*(K-1 < fN); % Everything ON
R2 = 2*(Son > K & (K-1)./Son > fN); % autonomous cells for Son > K
R3 = 3*(K./Son - 1 > fN); % Everything OFF
R4 = 4*(Son < K & K - Son < fN & (K-1)./Son > fN); % autonomous cells for Son < K

map = [1, 1, 1
    1, 0, 0
    1, 1, 0
    1, 0, 1
    0, 1, 0];
    
out = R1 + R2 + R3 + R4;
him = imagesc(K_vec, Son_vec, out);
tmp = map(1:max(max(out))+1,:);
colormap(tmp);
hold on
plot(K_vec, K_vec, 'k--')
hold off
set(gca,'ydir','normal')
colorbar;
