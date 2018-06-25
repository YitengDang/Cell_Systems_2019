%% Time evolution of lattice with visualization
clear all
close all
set(0,'defaulttextinterpreter', 'latex');
%%
% geometric parameters
Lx = 1;
n = 10; % nmax = L/R (square packing)
N = n^2; % total number of particles
%N = round(eta*L^2/(pi*R^2));
Ly = sqrt(3)/2*Lx;
r_av = sqrt(Lx*Ly/N)/2; % estimate for NND of random distribution (Clark & Evans, 1954)
rcell = 0.2;
R = rcell*Lx/(n+1); % disc radius
eta = N*pi*R^2/(Lx*Ly); %packing fraction 
if eta > pi/(2*sqrt(3)) % max. packing fraction
    disp('packing fraction too high! Abort evaluation');
    pause(10);
end

% Circuit Parameters
K = 16;
Con = 8;
noise = 0;
a0 = 0.5;
rcell = 0.2;

% (1) Markov MC
mcsteps = 10^4; %10^3;
[pos, dist, fN0] = initial_cells_random_markov_periodic(n, Lx, R, mcsteps);
% (2) random placement
%
%[pos, dist] = initial_cells_random_periodic_alt(n, Lx, Ly, R); % random config
%}

%% Initial configuration
% initial config
p = 0.7;
iniON = round(N*p);
t = 0;

hin = figure(1);
plot_handle = reset_cell_figure(hin, pos, rcell);
cells = zeros(N, 1);
cells(randperm(N, iniON)) = 1;
update_figure_periodic_scatter(plot_handle, cells, t);

%set(gca, 'XTick', [], 'YTick', []);
cells_hist = {};
cells_hist{end+1} = cells;

% Time evolution
changed = 1;
while changed
    % calculate new state
    pause(0.8);
    [cells_out, changed] = ...
        update_cells_noise(cells, dist, Con, K, a0, rcell*a0, noise);
    
    % update cells & figure
    t = t+1;
    cells = cells_out;
    cells_hist{end+1} = cells;
    update_figure_periodic_scatter(plot_handle, cells, t);
end

%% Nearest neighbour distances
%{
dist2 = reshape(dist(dist>0), [N-1 N]);
nnd = min(dist2);
nnd_av = mean(nnd);

% Histogram of NND
h=figure();
hold on
nbins = 20;
histogram(nnd, nbins, 'Normalization', 'pdf');
xlabel('NND');
ylabel('Probability');
title(sprintf('$$ \\langle NND \\rangle = %.2f$$, target $$= %.2f $$', nnd_av, r_av));
set(gca, 'FontSize', 24);
set(h, 'Position', [500 200 840 720]);

% target a0
plot([r_av r_av], [0 100], 'r--');

% Clark & Evans calculation
%{
rho = @(r) 2*pi*N/L^2*r.*exp(-pi*N/L^2*r.^2);
rvals = linspace(0, max(nnd)*1.2, 1000);
plot(rvals, rho(rvals), 'LineWidth', 2);
%xlim([0 1]);
%}
%}
%% Problem: distances do not match with those of original calculation
%[pos2,ex,ey] = init_cellpos_hex(n,n);
%dist2 = dist_mat(pos2,n,n,ex,ey);