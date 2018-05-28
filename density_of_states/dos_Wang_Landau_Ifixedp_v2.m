% ----Implements the Wang-Landau algorithm for calculating density of states----
% Short description: http://stp.clarku.edu/simulations/ising/wanglandau/index.html
% Then applies the Metropolis algorithm to obtain a normalized distribution.

% Article: D. P. Landau, Shan-Ho Tsai, and M. Exler,
% "A new approach to Monte Carlo simulations in statistical physics: Wang-Landau sampling," 
% Am. J. Phys. 72, 1294–1302 (2004).

% Cellular automaton model: calculate P(I|p), probability of obtaining I
% given value of p.

% To do: Introduce multiple walkers, parallelize
%--------------------------------------------------------
clear all;
close all;
clc;

%--Parameters--
% model parameters
L = 11; % size of lattice
N = L^2;
a0 = 0.5;
p=12/N;
%Rcell = 0.2*a0;
%Son = 8;
%K = 16;

% algorithm parameters
nbins = 30; % number of histogram bins
c1 = 1; f0 = exp(c1); % initial modification factor 
c2 = 4; ffinal = exp(10^(-c2)); % final modification factor
p_flat = 0.8; % flatness parameter
qsave = 0;
f = f0; % modification factor
%%
% ----Algorithm----
[dist, ~] = init_dist_hex(L, L);
iniON = round(p*N);
%{
% Initial random sampling to determine bin boundaries
cells = zeros(N,1); %random lattice
cells(randperm(N, iniON)) = 1;
[I, theta] = moranI(cells, a0*dist); % I of lattice
ntrials=10000;
allI = zeros(ntrials, 1);
for i=1:ntrials 
    cells = zeros(L^2,1);
    cells(randperm(N, iniON))=1;
    allI(i) = moranI(cells, a0*dist); % energy of lattice
end
%}
%epsilon = 0.02;
%edges=linspace(floor(min(allI)*1000)/1000, ceil(max(allI)*1000)/1000 + epsilon, nbins+1);

% Fixed bins for I
edges = linspace(-0.06, 0.14, nbins+1);
%%
% main loop
cells = zeros(N,1); %random lattice
cells(randperm(N, iniON)) = 1;
lg = dos_Wang_Landau_loop(cells, f, edges, ffinal, p_flat, nbins, dist, a0);

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

% Normalise g(E) by reshifting
lgn2 = lg - min(lg);

% Plot normalised g(E)
h1a=figure();
plot(bincenters, exp(lgn)/2^N, 'x-');
xlabel('I');
ylabel('P(I|p)');

% Plot g(E) normalised by reshifting
h1b=figure();
plot(bincenters, exp(lgn2)/sum(exp(lgn2)), 'x-');
xlabel('I');
ylabel('P(I|p)');
%%
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
qsave = 0;
if qsave
    % save figure
    % h1a: exact normalisation, h1b: normalisation by reshifting
    fname = strrep(sprintf('WL_norm_N%d_n0_%d_a0_%.2f_f0exp%d_ffin_e10e-%d_pflat%.1f_%dbins_fixed%.2f',...
        N, p*N, a0, c1, c2, p_flat, nbins), '.', 'p');
    out_fig = fullfile(pwd, 'figures', 'omega_pI', 'Wang-Landau probI_fixedp', ...
        strcat(fname,'.pdf')); % filename figure
    save_figure_pdf(h1b, 10, 6, out_fig);
    % export data
    out_file = fullfile(pwd, 'data', 'dos', 'probI_fixedp', strcat(fname, '.mat') ); % filename
    close all;
    clear dist allI; %clear large vars that are not needed
    save(out_file);
end
%% Calculate Omega(p,I)
% Omega(p,I) = P(I|p)Omega(p)
omegap = nchoosek(N,iniON);
omegapI = omegap/2^N*exp(lgn2);
h2=figure();
plot(bincenters, omegapI)
xlabel('$$I$$');
ylabel('$$\Omega(p,I)$$');
