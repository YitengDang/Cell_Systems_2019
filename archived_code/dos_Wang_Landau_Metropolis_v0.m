% ----Implements the Wang-Landau algorithm for calculating density of states----
% Short description: http://stp.clarku.edu/simulations/ising/wanglandau/index.html
% Then applies the Metropolis algorithm to obtain a normalized distribution.

% Article: D. P. Landau, Shan-Ho Tsai, and M. Exler,
% "A new approach to Monte Carlo simulations in statistical physics: Wang-Landau sampling," 
% Am. J. Phys. 72, 1294–1302 (2004).

% Cellular automaton model: calculate P(I | p), probability of obtaining I
% given value of p.

% CONCLUSION: DOES NOT WORK BECAUSE METROPOLIS SAMPLES OVER MOST LIKELY
% STATES AND MISSES OUT MANY RARE STATES
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
p = 0.5;
iniON = round(p*N);

% algorithm parameters
nbins = 30; % number of histogram bins
c1=1; f0 = exp(c1); % initial modification factor 
c2=2; ffinal = exp(10^(-c2)); % final modification factor
p_flat=0.5; % flatness parameter
qsave = 1;

%--Initialization--
% model
[dist, ~] = init_dist_hex(L, L);

% algorithm
edges = linspace(-0.1,0.3,nbins+1); % I value bin edges
f = f0; % modification factor 
%g = ones(1,nbins); % estimated g(E), updated in loop
lg = ones(1,nbins); % estimated log(g(E))
%energy 'histogram' initialised inside loop

cells = zeros(N,1); %random lattice
cells(randperm(N, iniON)) = 1;
[I, theta] = moranI(cells, dist); % I of lattice
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
        % flip two random spins: 1 ON and 1 OFF
        cellsp = cells;
        ONcells = find(cells == 1); x1=randi(length(ONcells)); cellsp(ONcells(x1))=0; % turn ON cell OFF
        OFFcells = find(cells == 0); x2=randi(length(OFFcells)); cellsp(OFFcells(x2))=1; % turn OFF cell ON
        Ip = moranI(cellsp, dist);
        bin = dos_Wang_Landau_binning(I, edges); % bin of original system
        binp = dos_Wang_Landau_binning(Ip, edges); % bin of flipped system
        %disp(h); disp(hp);
        %disp(g(E)); disp(g(Ep));
        %disp(lg(E)); disp(lg(Ep));
        prob = min(exp(lg(bin)-lg(binp)), 1); %acceptance probability min(g()/g(E2), 1)
        %p = min(g(E)/g(Ep), 1); %acceptance probability min(g()/g(E2), 1)
        if rand() < prob
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
set(0,'defaulttextinterpreter','latex');

h=figure();
bincenters = (edges(1:end-1)+edges(2:end))/2;
plot(bincenters, lg, 'x-');
xlabel('I'); %, 'Interpreter', 'latex');
%ylabel('$$P(I|p)$$', 'Interpreter', 'latex');
ylabel('$$\log{P(I|p)}$$'); %, 'Interpreter', 'latex');

% Normalise g(E)
 lgn = lg - log(sum(exp(lg))) + N*log(2); %Method 1 (may lead to Inf)
% lgn = lg - lg(1) + log(2); % Method 2, works only if first bin contains only the ground state

% Plot normalised g(E)
h1=figure();
plot(bincenters, exp(lgn)/2^N, 'x-');
xlabel('I');
ylabel('P(I|p)');

% Save figure
%{
if save
    fname = strrep(sprintf('Wang_Landau_N%d_p0_%.2f_f0exp%d_ffinalexp10e-%d_pflat_%.1f_%dbins',...
        N, p, c1, c2, p_flat, nbins), '.', 'p');
    out_file = fullfile(pwd, 'figures', 'omega_pI', 'Wang-Landau', ...
        strcat(fname,'.pdf')); % filename
    save_figure_pdf(h, 10, 6, out_file);    
end
%}
%%
if qsave
    fname = strrep(sprintf('Wang_Landau_normalised_N%d_p0_%.2f_f0exp%d_ffinalexp10e-%d_pflat%.1f_%dbins',...
        N, p, c1, c2, p_flat, nbins), '.', 'p');
    out_fig = fullfile(pwd, 'figures', 'omega_pI', 'Wang-Landau', ...
        strcat(fname,'.pdf')); % filename figure
    %save_figure_pdf(h1, 10, 6, out_fig);
    % export data
    out_file = fullfile(pwd, 'data', 'dos', 'probI_fixedp', strcat(fname, '.mat') ); % filename
    close all;
    clear dist; %clear large vars that are not needed
    save(out_file);
end
%%
% After flat distribution is obtained, sample over this distribution using
% Metropolis (almost the same, but acceptance probability reversed!)

nruns = 100000; %datapoints to take
P = zeros(1,nbins);
for i=1:nruns 
        % flip two random spins: 1 ON and 1 OFF
        cellsp = cells;
        ONcells = find(cells == 1); x1=randi(length(ONcells)); cellsp(ONcells(x1))=0; % turn ON cell OFF
        OFFcells = find(cells == 0); x2=randi(length(OFFcells)); cellsp(OFFcells(x2))=1; % turn OFF cell ON
        Ip = moranI(cellsp, dist);
        bin = dos_Wang_Landau_binning(I, edges); % bin of original system
        binp = dos_Wang_Landau_binning(Ip, edges); % bin of flipped system
        %disp(h); disp(hp);
        %disp(g(E)); disp(g(Ep));
        %disp(lg(E)); disp(lg(Ep));
        prob = min(exp(lg(binp)-lg(bin)), 1); %acceptance probability min(g()/g(E2), 1)
        if rand() < prob
            cells = cellsp;
            I = Ip;
            bin = binp;
        end
        %disp(lattice);
        %disp(strcat('new h=', num2str(h)));
        
        %bin energies
        P(bin) = P(bin) + 1;
end
%% Plot final P(I)
h=figure();
bincenters = (edges(1:end-1)+edges(2:end))/2;
plot(bincenters, P/nruns, 'x-');
xlabel('I', 'Interpreter', 'latex');
ylabel('$$P(I|p)$$', 'Interpreter', 'latex');
