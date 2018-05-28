% ----Implements the Wang-Landau algorithm for calculating density of states----
% Short description: http://stp.clarku.edu/simulations/ising/wanglandau/index.html
% 
% Article: D. P. Landau, Shan-Ho Tsai, and M. Exler,
% "A new approach to Monte Carlo simulations in statistical physics: Wang-Landau sampling," 
% Am. J. Phys. 72, 1294–1302 (2004).

% Cellular automaton model: calculate P(I | p), probability of obtaining I
% given value of p.
%--------------------------------------------------------
clear all;
close all;
clc;

%--Parameters--
% model parameters
L = 11; % size of lattice
N = L^2;
%a0 = 0.5;
%Rcell = 0.2*a0;
%Son = 8;
%K = 16;
p = 0.3;
iniON = round(p*N);

% algorithm parameters
nbins = 30; % number of histogram bins
c1=1; f0 = exp(c1); % initial modification factor 
c2=2; ffinal = exp(10^(-c2)); % final modification factor
p_flat=0.5; % flatness parameter
save = 1;

%--Initialization--
% model
[dist, ~] = init_dist_hex(L, L);

% algorithm
edges = linspace(-0.15,0.5,nbins+1); % energy bin edges
f = f0; % modification factor 
%g = ones(1,nbins); % estimated g(E), updated in loop
lg = ones(1,nbins); % estimated log(g(E))
%H = ones(1,nbins); % energy 'histogram', initialised inside loop

cells = zeros(N,1); %random lattice
cells(randperm(N, iniON)) = 1;
I = moranI(cells, dist); % I of lattice
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
        % choose position of two random spins: 1 ON and 1 OFF
        %x = randi(N, 1, 2);
        cellsp = cells;
        %cellsp(x) = ~cellsp(x); %trial flip
        ONcells = find(cells == 1); x1=randi(length(ONcells)); cellsp(ONcells(x1))=0; % turn ON cell OFF
        OFFcells = find(cells == 0); x2=randi(length(OFFcells)); cellsp(OFFcells(x2))=1; % turn OFF cell ON
        Ip = moranI(cellsp, dist);
        bin = Wang_Landau_binning(I, edges); % bin of original system
        binp = Wang_Landau_binning(Ip, edges); % bin of flipped system
        %disp(h); disp(hp);
        %disp(g(E)); disp(g(Ep));
        %disp(lg(E)); disp(lg(Ep));
        p = min(exp(lg(bin)-lg(binp)), 1); %acceptance probability min(g()/g(E2), 1)
        %p = min(g(E)/g(Ep), 1); %acceptance probability min(g()/g(E2), 1)
        if rand() < p
            cells = cellsp;
            I = Ip;
            bin = binp;
        end
    %disp(lattice);
    %disp(strcat('new h=', num2str(h)));
    
    %-update ge and energies
    %g(E) = f*g(E);
    lg(bin) = lg(bin) + log(f);
    H(bin) = H(bin) + 1;
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
% lgn = lg - lg(1) + log(2); % Method 2, works only if first bin contains only the ground state

h=figure();
bincenters = (edges(1:end-1)+edges(2:end))/2;
plot(bincenters, lg, 'x-');
xlabel('I', 'Interpreter', 'latex');
%ylabel('$$P(I|p)$$', 'Interpreter', 'latex');
ylabel('$$\log{P(I|p)}$$', 'Interpreter', 'latex');

% Save figure
%save=1;
if save
    fname = strrep(sprintf('Wang_Landau_N%d_p0_%.2f_f0_exp%d_ffinal_exp10emin%d_pflat_%.1f_%dbins_edges%.2fto%.2f',...
        N, p, c1, c2, p_flat, nbins, edges(1), edges(end)), '.', 'p');
    out_file = fullfile(pwd, 'figures', 'omega_pI', 'probI_fixedp', ...
        strcat(fname,'.pdf')); % filename
    save_figure_pdf(h, 10, 6, out_file);
end
%% --------------Random lattice trials to fix system parameters------------------
% Average sampling time given by coupon problem

%--Parameters--
% model parameters
L = 11; % size of lattice
N = L^2;
%a0 = 0.5;
%Rcell = 0.2*a0;
%Son = 8;
%K = 16;
p = 0.5;
iniON = round(p*N);

%--Initialization--
[dist, ~] = init_dist_hex(L, L);

% Plot histogram I
ntrials = 100000;
Ivals = zeros(1,ntrials);
for i=1:ntrials
    %%cells = zeros(N,1); %random lattice
    %%cells(randperm(N, iniON)) = 1;
    cells = randi(2,N,1)-1;
    Ivals(i) = moranI(cells, dist); % I of lattice
    %disp(lattice);
end
h2=figure();
histogram(Ivals, 'Normalization', 'probability');
xlabel('I');
ylabel('probability');
%}

% Save figure
save=1;
if save
    fname = sprintf('N%d_n0%d_p0%.2f_%dtrials', N, iniON, p, ntrials);
    out_file = fullfile(pwd, 'figures', 'omega_pI', 'probI_fixedp_estimate', ...
        strcat(fname,'_random_sampling.pdf')); % filename
    save_figure_pdf(h2, 10, 6, out_file);
end
%% Try to find lattice configurations with desired energies
h_target = 12:2:22;
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