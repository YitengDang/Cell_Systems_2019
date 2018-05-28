% ----Implements the Wang-Landau algorithm for calculating density of states----
% Short description: http://stp.clarku.edu/simulations/ising/wanglandau/index.html
% 
% Article: D. P. Landau, Shan-Ho Tsai, and M. Exler,
% "A new approach to Monte Carlo simulations in statistical physics: Wang-Landau sampling," 
% Am. J. Phys. 72, 1294–1302 (2004).

% Calculates P(H) for the multicellular Hamiltonian H
%--------------------------------------------------------
clear all;
close all;
clc;

%--Parameters--
% model parameters
L = 3; % size of lattice
a0 = 0.5;
Rcell = 0.2*a0;
Son = 8;
K = 16;

% algorithm parameters
nbins = 10; % number of histogram bins
f0 = exp(2); % initial modification factor 
ffinal = exp(10^(-2)); % final modification factor
p_flat=0.4; % flatness parameter
save = 1;

%--Initialization--
% model
[dist, pos] = init_dist_hex(L, L);

% algorithm
edges = linspace(12,22,nbins+1); % energy bin edges
f = f0; % modification factor 
%g = ones(1,nbins); % estimated g(E), updated in loop
lg = ones(1,nbins); % estimated log(g(E))
%H = ones(1,nbins); % energy 'histogram', initialised inside loop

cells = randi(2,L^2,1)-1; %random lattice of cells (stored as linear array)
h = hamiltonian(cells, dist, Son, K, a0, Rcell); % energy of lattice
%disp(lattice);
%disp(h);
% function binning(E, bins) determines which bin the energy E falls in
H = zeros(1,nbins); 
%%
% ----Algorithm----
% Introduce multiple walkers, parallelize

% main loop
iteration=1;

while f > ffinal
disp(strcat('Iteration:', num2str(iteration)));
H = zeros(1,nbins); % initialize/reset histogram 
% --inner loop
%for i=1:1000
flat=0;

while ~flat
    %disp('check');
    for i=1:10000 % Do 10000 MC steps before checking
        x = randi(L^2); % choose position of random spin
        cellsp = cells;
        cellsp(x) = ~cellsp(x); %trial flip
        hp = hamiltonian(cellsp, dist, Son, K, a0, Rcell);
        E = dos_Wang_Landau_binning(h, edges); % energy bin of original system
        Ep = dos_Wang_Landau_binning(hp, edges); % energy bin of flipped system
        %disp(h); disp(hp);
        %disp(g(E)); disp(g(Ep));
        %disp(lg(E)); disp(lg(Ep));
        p = min(exp(lg(E)-lg(Ep)), 1); %acceptance probability min(g()/g(E2), 1)
        %p = min(g(E)/g(Ep), 1); %acceptance probability min(g()/g(E2), 1)
        if rand() < p
            cells = cellsp;
            h = hp;
            E = Ep;
        end
    %disp(lattice);
    %disp(strcat('new h=', num2str(h)));
    
    %-update ge and energies
    %g(E) = f*g(E);
    lg(E) = lg(E) + log(f);
    H(E) = H(E) + 1;
    %disp(g);
    %disp(H);
    end
    
    %-check flatness criterion (all H(E) > p*<H(E)>)
    flat = all(H > mean(H)*p_flat); %flatness criterion
    %disp(flat);
    %end
    % --end inner loop
end
f = f^(1/2);
disp(strcat('f=',num2str(f)));
iteration = iteration+1;
end
%% Plot final g(E)
% Normalise g(E)
% lgn = lg - log(sum(exp(lg))) + N*log(2); %Method 1 (may lead to Inf)
lgn = lg - lg(1) + log(2); % Method 2, works only if first bin contains only the ground state

h=figure();
bincenters = (edges(1:end-1)+edges(2:end))/2;
plot(bincenters, lgn, 'x-');
xlabel('Energy', 'Interpreter', 'latex');
ylabel('$$\log{g(E)}$$', 'Interpreter', 'latex');
%plot(bincenters, log(g), 'x-');

% Save figure
fname = sprintf('Wang_Landau_L%d_f_exp1_ffinal_exp10tomin8_pflat_%d_normalised', L, round(10*p_flat));
if save
    out_file = fullfile(pwd, 'data', 'dos', fname);
    save_figure_pdf(h, 10, 6, out_file);
end
%%
%figure();
%bar(bincenters,lg,1);
%xlabel('Energy');
%ylabel('$$\log{g(E)}$$', 'Interpreter', 'latex');

%% Random lattice trials to fix system parameters
% Average sampling time given by coupon problem
ntrials = 10000;
hvals = zeros(1,ntrials);
for i=1:ntrials
    cells = randi(2,L^2,1)-1; %random lattice
    hvals(i) = hamiltonian(cells, dist, Son, K, a0, Rcell); % energy of lattice
    %disp(lattice);
end
figure();
histogram(hvals);
%}

%% Try to find lattice configurations with desired energies
h_target = [12:2:22];
margin = 1;
target_cells = zeros(length(h_target), L^2);

for i=1:length(h_target)
    found = 0;
    while ~found
        cells = randi(2,L^2,1)-1; %random lattice
        h_rand = hamiltonian(cells, dist, Son, K, a0, Rcell); % energy of lattice
        if abs(h_rand-h_target(i)) < margin
            target_cells(i,:) = cells;
            disp(h_target(i));
            disp(cells);
            found = 1;
        end
    end
end