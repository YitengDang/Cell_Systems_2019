% Calculates the average value of the multicellular hamiltonian in K, C_ON
% space for fixed values of a0
clear variables
close all
clc
%%
% Fixed variables
gridsize = 11;
N = gridsize^2;
a0 = 1.5;
Rcell = 0.2*a0;

% calculate fN
[dist, pos] = init_dist_hex(gridsize, gridsize);
dist_vec = a0*dist(1,:);
r = dist_vec(dist_vec>0); % exclude self influence
fN = sum(sinh(Rcell)*sum(exp(Rcell-r)./r)); % calculate signaling strength

% Load Omega(p,I) data
prange = (0:floor(N/2))/N;
Irange = linspace(0.06, 0.14, 20);
% NOT IMPLEMENTED YET!
%%
% Compute <H>
K = linspace(1,20,100);
Son = linspace(1,30,100);
h = @(p, I, K, Son) -0.5*(Son-1).*(1 + 4*fN.*p.*(1-p).*I + fN*(2*p-1).^2) ...
                -(2*p-1).*(0.5*(Son+1)*(1+fN) - K);
[pvals, Ivals] = meshgrid(prange,Irange);
omegapI = ones(size(pvals));
havg = zeros(100,100);
for i=1:100
    for j=1:100
        thisK = K(i);
        thisSon = Son(j);
        havg(i,j) = sum(sum(h(pvals, Ivals, thisK, thisSon).*omegapI));
    end
end
%%
% Plot <h>(K, Con)
set(0, 'defaulttextinterpreter', 'latex');
figure();
imagesc(K, Son, havg);
c = colorbar;
set(gca,'FontSize', 20);
xlabel('$$K$$', 'FontSize', 24)
ylabel('$$C_{ON}$$', 'FontSize', 24)
ylabel(c, 'h=H/N', 'FontSize', 24)
%% Compute estimate for W 
% for tuning K

% protocol 4
tsteps = 10;
Kprot = linspace(12, 2, tsteps+1);
Cprot = 2*Kprot/(1+fN) - 1;
a0prot = 1.5;
dK = Kprot(2:end)-Kprot(1:end-1);

% define weights
pvals = (0:N)/N;
% Calculate binomial coefficients (possibly not exact)
% first check if file exists already
fileid = sprintf('nchoosek_n%d', N);
fname = fullfile(pwd,'data','nchoosek', strcat(fileid, '.mat'));
if exist(fname)==2
    disp('exists!');
    load(fname);
else
    disp('Doesnt exist!');
    binom = zeros(1,N+1);
    for j=1:N+1  
        disp(j);
        % if not, calculate
        binom(j) = nchoosek(N,j-1);
    end
    save(fname, 'binom');
end

% calculate fN
dist_vec = a0prot*dist(1,:);
Rcell = 0.2*a0prot;
fN = sum(sinh(Rcell)*sum(exp(Rcell-r)./r)); % calculate signaling strength

% Define forces
syms p I
fK = @(p) -(2*p - 1);
%%
% define work
wK = sum(dK)*(sum((2*pvals).*(binom/2^N)) - 1);
%wCon = sum(sum(-fCon(pvals, Ivals).*dCon)*pp(pvals));

%%
fCon = @(p,I) 0.5*(1 + fN*p*(1-p)*I - (2*p-1)^2*fN) + (2*p-1)*(1+fN)/2;
