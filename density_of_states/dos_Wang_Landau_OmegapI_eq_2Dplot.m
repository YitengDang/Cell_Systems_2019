% Plots results of Wang-Laudau algorithm for Omega(p,I) in 2D plots
clear all;
close all;
clc;
%% Load files
L = 11; % size of lattice
N = L^2;
a0 = 0.5;
K=16;
Son=8;
Rcell=0.2*a0;
%
n = 11:61;
p=n/N;
nbins = 30; % number of histogram bins
c1=1; f0 = exp(c1); % initial modification factor 
c2=4; ffinal = exp(10^(-c2)); % final modification factor
p_flat=0.8; % flatness parameter
qsave = 1;

g_all = zeros(nbins, numel(n));
for i=1:numel(n)
    this_n = n(i);
    fileid = strrep(sprintf('WL_norm_N%d_n0_%d_a0_%.2f_f0exp%d_ffin_e10e-%d_pflat%.1f_%dbins_fixed%.2f',...
        N, this_n, a0, c1, c2, p_flat, nbins), '.', 'p');
    path = fullfile(pwd, 'data', 'dos', 'probI_fixedp', strcat(fileid, '.mat') ); % filename
    load(path, 'lgn2', 'edges');
    g_all(:, i) = exp(lgn2)/sum(exp(lgn2));
end
%% --Plot Omega_E(p,I)--
% Calculate data for Omega_E(p)
% Initiate lattice
[dist, pos] = init_dist_hex(L, L);

% Calculate fN and gN
dist_vec = a0*dist(1,:);
r = dist_vec(dist_vec>0); % exclude self influence
fN = sinh(Rcell)*sum(exp(Rcell-r)./r); % calculate signaling strength
gN = sinh(Rcell)^2*sum(exp(2*(Rcell-r))./(r.^2)); % calculate signaling strength

% Calculate Ponon and Poffoff
mup = ((Son-1).*p + 1)*fN;
sigmap = (Son-1)^2.*p.*(1-p)*gN; 
Ponon = 1-normcdf(K-Son, mup, sigmap);
Poffoff = normcdf(K-1, mup, sigmap);

% Calculate binomial coefficients (possibly not exact)
% first check if file exists already
fileid = sprintf('%d', N);
fname = fullfile(pwd,'data','nchoosek', strcat(fileid, '.mat'));
if exist(fname)==2
    disp('exists!');
    load(fname);
else
    disp('Doesnt exist!');
    nck = zeros(N+1,1);
    for j=1:N+1  
        disp(j);
        % if not, calculate
        nck(j) = nchoosek(N,j-1);
    end
    save(fname, 'binom');
end
binom = nck(n);

% Calculate steady state pattern
omegap = Ponon.^(p*N).*Poffoff.^(N*(1-p)).*binom';
figure();
plot(p,omegap)
%%
% Plot figures
set(0, 'defaulttextinterpreter', 'latex');
Icenters = (edges(1:end-1)+edges(2:end))/2;
h1=figure();
omegapI = g_all.*repmat(omegap, nbins, 1);
colormap('hot');
X=sort([n N-n])/N;
imagesc(X, Icenters, log([omegapI fliplr(omegapI)]) );
c=colorbar;
c.Label.String = 'log(\Omega(p,I))';
set(gca,'FontSize', 14);
set(gca,'YDir','normal');
xticks(0.1:0.1:0.9);
yticks(edges(1:3:end));
xlabel('p');
ylabel('I');

%% Save heat map
if qsave
    fname = strrep(sprintf('OmegapI_WL_norm_N%d_a0_%.2f_f0exp%d_ffin_e10e-%d_pflat%.1f_%dbins_fixed%.2f',...
        N, a0, c1, c2, p_flat, nbins), '.', 'p');
    out_fig = fullfile(pwd, 'figures', 'dos', 'Wang-Landau probI_fixedp', ...
        strcat(fname,'.pdf')); % filename figure
    save_figure_pdf(h1, 10, 6, out_fig);
end