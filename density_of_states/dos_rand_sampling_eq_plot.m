% Calculates the density of equilibrium states P_E(p,I) from P_E(p,I)=P(I|p)P_E(p)
% P_E(p) is known analytically (approximation)
% P(I|p) is determined through simulations (requires
% density_of_states_random.m)
% Assume same distribution for equilibrium states, Omega_E(p,I) =
% P(I|p)Omega_E(p)
% v2: P(I|p) from random sampling now for different fixed p

clear all
close all
L=11;
N=L^2;
p=(0:N)/N; 
a0 = 1.5;
Rcell = 0.2*a0;
K = 15;
Son = 20;
qsave = 0;
%% ---(1) Plot Omega(p)---
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

% Calculate steady state pattern
omegap = Ponon.^(p*N).*Poffoff.^(N*(1-p)).*nck';

% Plot Omega(p)
set(0,'defaulttextinterpreter','latex');
figure();
plot(p, omegap, 'x-');
xlabel('p');
ylabel('$$\Omega_E(p)$$');%, 'Interpreter', 'latex');
%semilogy(p, omegap);
%% ---(2) Plot Omega(p, I)---
% ---P(I|p) from random simulations, P(p) analytical---
% (a) Import data from calculating P(p,I) 
ntrials = 10000;
nbins = 100;
epsilon = 0.01;

fileid = strrep(sprintf('OmegapI_rand_fixedp_N%d_runs%d_bins%d_margin%.2f',...
        L^2, ntrials, nbins, epsilon), '.','p');
fname = fullfile(pwd, 'data', 'omega_pI', 'random_sampling', ...
        strcat(fileid,'.mat'));
load(fname);    

% Plot Omega(I) histogram from loaded data
%{
figure();
nbins = 100; %#bins for binning data AND for histogram
epsilon = 0.02;
Iedges = linspace(min(Ivals)-epsilon, max(Ivals)+epsilon, nbins+1);
Ibincenters = (Iedges(1:end-1)+Iedges(2:end))/2;
plot(Ibincenters, Ibincounts/ntrials, 'x-');
%}
% (b) Calculate Omega(p,I)
Icondp = Ibincounts/ntrials; % Conditional probability P(I|p), rows=p, columns=I 

% multiply each column of omegap with a row of Icondp to get P(I|p) P(p)
% for that value of p
omegapI = zeros(N+1, nbins);
for i=1:N+1
    omegapI(i,:) = Icondp(i,:)*omegap(i);
end

% (b2) set values below 1 equal to 0
idx = find(omegapI < 1);
omegapI(idx) = 0;
%%
% (c) Plot as surface
% mesh parameters
%nbins = 100; %#bins for binning data AND for histogram
%epsilon = 0.01;
%Iedges = linspace(min(Ivals)-epsilon, max(Ivals)+epsilon, nbins+1);
Icenters = (edges(1:end-1)+edges(2:end))/2;
[X,Y] = meshgrid(p, Icenters); 

% plot surface
%set(0,'defaulttextinterpreter','latex');
h1=figure();
surf(X,Y,log(omegapI'));
%mesh(log(OmegapI));
xlabel('p');
ylabel('I');
zlabel('$$\log{\Omega_E(p,I)}$$'); %, 'Interpreter', 'latex');
%%
% (d) Plot as heat map
h2=figure();
colormap('summer');
imagesc(p, Icenters, log(omegapI'));
c=colorbar;
c.Label.String = 'log(Omega_E(p,I))';
set(gca,'FontSize', 14);
set(gca,'YDir','normal')
xlabel('p');
ylabel('I');

%% Compare with Hamiltonian
E = @(p, I) -0.5*(Son-1)*(4*fN.*p.*(1-p).*I + fN*(2*p-1).^2) ...
                -(2*p-1).*(0.5*(Son+1)*fN - K);
            
% Heat map
[pmesh,Imesh] = meshgrid(p, Icenters);
Hvals = E(pmesh,Imesh);

%{
% Heat map
figure();
colormap('hot');
imagesc(p, Icenters, Hvals);
c=colorbar;
c.Label.String = 'H(p,I)';
set(gca,'FontSize', 14);
xlabel('p');
ylabel('I');
%}
% Contour map
h3=figure();
contourf(p, Icenters, Hvals,'LineStyle', 'none')
colormap('summer')
c = colorbar;
set(gca,'FontSize', 14);
xlabel('p');
ylabel('I');

%%
% Save figures
if qsave
    fileid = strrep(sprintf('dos_eq_states_Isim_N%d_a0_%.1f_Con%d_K%d_ntrials%d_test', N, a0, Son, K, ntrials), '.', 'p');
    %out_name1 = fullfile(pwd, 'figures', 'omega_pI', 'p_exact_I_sim', strcat(fileid, '_surf_map.pdf'));
    %save_figure_pdf(h1, 10, 6, out_name1);
    out_name2 = fullfile(pwd, 'figures', 'omega_pI', 'p_exact_I_sim', strcat(fileid, '_heat_map.pdf'));
    save_figure_pdf(h2, 10, 6, out_name2);
    out_name3 = fullfile(pwd, 'figures', 'omega_pI', 'p_exact_I_sim', strcat(fileid, '_hamiltonian.pdf'));
    save_figure_pdf(h3, 10, 6, out_name3);
end