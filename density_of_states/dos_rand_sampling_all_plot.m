% Calculates the density of states P(p,I) from P(p,I)=P(I|p)P(p)
% P_E(p) is known analytically (approximation)
% P(I|p) is determined through simulations (requires
% dos_wang_landau_Ifixedp_v2.m)
% Assume same distribution for equilibrium states, Omega_E(p,I) =
% P(I|p)Omega_E(p)
% v2: P(I|p) from random sampling now for different fixed p
clear all
close all
clc
d1 = digits(64);

L = 11;
N = L^2;
p = (0:N)/N;
qsave = 0;
%% ---(1) Plot Omega(p)---
% Initiate lattice
[dist, pos] = init_dist_hex(L, L);

% Calculate binomial coefficients (possibly not exact)
% first check if file exists already
fileid = sprintf('nchoosek_n%d', N);
fname = fullfile(pwd,'data','nchoosek', strcat(fileid, '.mat'));
if exist(fname)==2
    %disp('exists!');
    load(fname);
else
    [out, msg] = mkdir('data','nchoosek');
    %disp('Doesnt exist!');
    binom = zeros(1,N+1);
    for j=1:N+1  
        disp(j);
        % if not, calculate
        binom(j) = nchoosek(N,j-1);
    end
    save(fname, 'binom');
end

% Omega(p) 
omegap = binom;

% Plot Omega(p)
set(0,'defaulttextinterpreter','tex'); %works for all future plots
figure();
plot(p, omegap, 'x-');
xlabel('p');
ylabel('$$\Omega(p)$$'); %, 'Interpreter', 'latex');
%semilogy(p, omegap);
%% ---(2) Plot Omega(p, I)---
% ---P(I|p) from random simulations, P(p) analytical---
% (a) Import data from calculating P(p,I) 
ntrials = 10000;
nbins = 100;
epsilon = 0.01;

%{
fileid = strrep(sprintf('probIfixedp_rand_fixedp_N%d_runs%d_bins%d_margin%.2f',...
        L^2, ntrials, nbins, epsilon), '.','p');
fname = fullfile('H:\My Documents\Multicellular automaton', 'data', 'dos', 'probI_fixedp', ...
        strcat(fileid,'.mat'));
%}
fname = 'H:\My Documents\Multicellular automaton\figures\dos\random_sampling\probIrandp_rand_N225_runs1000000_bins100';
load(fname);    
%Ibincounts = numel(Iedges)-1;

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
Iedges = linspace(min(Ivals)-epsilon, max(Ivals)+epsilon, nbins+1);
Ibincenters = (Iedges(1:end-1)+Iedges(2:end))/2;
[X,Y] = meshgrid(p, Ibincenters);

% plot surface
%set(0,'defaulttextinterpreter','latex');
h1=figure();
surf(X,Y,log(omegapI'));
%mesh(log(OmegapI));
xlabel('p');
ylabel('I');
zlabel('$$\log{\Omega(p,I)}$$'); %, 'Interpreter', 'latex');
%%
% (d) Plot as heat map
h2=figure();
colormap('hot');
imagesc(p, Ibincenters, log(omegapI'));
c=colorbar;
c.Label.String = 'log(Omega(p,I))';
set(gca,'FontSize', 14);
xlabel('p');
ylabel('I');
set(gca,'YDir','normal')

%%
% Save figures
if qsave
    fileid = strrep(sprintf('dos_all_states_Isim_N%d_ntrials%d', N, ntrials), '.', 'p');
    out_name1 = fullfile(pwd, 'figures', 'dos', 'p_exact_I_sim', strcat(fileid, '_surf_map_p_fixed.pdf'));
    save_figure_pdf(h1, 10, 6, out_name1);
    out_name2 = fullfile(pwd, 'figures', 'dos', 'p_exact_I_sim', strcat(fileid, '_heat_map_p_fixed.pdf'));
    save_figure_pdf(h2, 10, 6, out_name2);
end