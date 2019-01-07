%% Time evolution of lattice with visualization
clear all
close all
clc

% geometric parameters
L = 2;
R = 0.02; % disc radius
n = 12; % nmax = L/R (square packing)
N = n^2; % total number of particles
%N = round(eta*L^2/(pi*R^2));
r_av = sqrt(L^2/N)/2; % estimate for NND of random distribution (Clark & Evans, 1954)
eta = N*pi*R^2/L^2; %packing fraction 
if eta > pi/(2*sqrt(3)) % max. packing fraction
    disp('packing fraction too high! Abort evaluation');
    pause(10);
end

% circuit parameters
K = 16; 
S = 8;

% Initial configuration
cells = ones(N, 1);
[pos, dist, rej] = initial_cells_random(N, L, R);
fprintf('Final rejections: %d \n', rej);

% draw figure
hin = figure(1);
update_figure(hin, pos, N, L, R, cells)

% time steps
tmax = 10; 

% Evolution
cells_hist = {};
cells_hist{end+1} = cells;

for i=1:tmax
    pause(1);
    [cells_out, changed] = update_cells(cells, dist, S, K, R);
    cells = cells_out;
    cells_hist{end+1} = cells;
    update_figure(hin, pos, N, L, R, cells);
    if sum(cells)==0
        break %if all cells are dead, stop
    end
end
%% Nearest neighbour distances
dist2 = reshape(dist(dist>0), [N-1 N]);
nnd = min(dist2);
nnd_av = mean(nnd);

h=figure();
hold on
histogram(nnd, 'Normalization', 'pdf');
xlabel('NND');
ylabel('Probability');
title(sprintf('%s NND %s = %.2f', '\langle', '\rangle', nnd_av));
set(gca, 'FontSize', 24);
set(h, 'Position', [500 200 840 720]);

% Clark & Evans calculation
rho = @(r) 2*pi*N/L^2*r.*exp(-pi*N/L^2*r.^2);
rvals = linspace(0, max(nnd)*1.2, 1000);
plot(rvals, rho(rvals), 'LineWidth', 2);
%xlim([0 1]);