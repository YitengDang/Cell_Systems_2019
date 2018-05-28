% Implements the Gillespie (stochastic stimulation) algorithm for diffusion
% on a hexagonal lattice.
 
% 1. Draw random number r
% 2. Propensity function alpha0 = k+gamma*X + 6*d
% 3. Reaction time 1/alpha0*log(r)
% 4. Adjust new concentration based on value of r.
    % -Choose values for parameters
%%
clear variables
close all
rng('shuffle');

% parameters
% k/gamma ~ number of molecules
% Aim for d*gamma ~ 1 to get lambda=sqrt(D/gamma) ~ 1.
L = 11; %gridsize
d = 10^3; % diffusion rate
k = 10^3; % production rate
gamma = 10^(-3); % degradation rate

% initialisation
X0 = 0;
X = X0*ones(L, L); % initial background concentration on lattice
cell = [1 1]; % secreting cell
X(1,1) = X0; % start with higher concentration on cell
t = []; % times
nsteps = 10^3; % simulation steps

% Define functions for moving in positive and negative directions
syms n
plus = @(n) n+1 - (n==L)*L;
minus = @(n) n-1 + (n==1)*L;
%%
% ---Keep track of choices (for troubleshooting)
events = zeros(1,8);
choice = zeros(L, L);
% ----------------------------------------------


for i=1:nsteps
% 1. Draw random number r
r = rand();
% 2. Propensity function alpha0 = k+gamma*X + 6*d
%alpha0 = 6*d;
alpha0 = k+gamma*sum(sum(X))+6*d;
% 3. Reaction time 1/alpha0*log(r)
t(end+1) = -1/alpha0*log(r);

% 4. Adjust new concentration based on value of r.
% (a) find cell that is not empty X(x,y) ~= 0, restore lattice if necessary
if sum(sum(X)) == 0  % if the lattice is empty
    disp('[X] = 0 on lattice!');
    X(1,1) = 1; %restore lattice by production
    t(end+1) = -1/k*log(rand()); % 
    x = 1; y = 1;
else %pick random particle of non-empty site
    r2 = rand();
    prob = X/sum(sum(X)); % probabilities proportional to site occupancy
    probc = cumsum([0, reshape(prob, 1, L^2)]); % cumulative probabilities
    idx = sum(r2 > probc); % linear index of chosen cell
    [x,y] = ind2sub([L,L], idx);
end

if r <= 6*d/alpha0 % diffusion
    %disp('opt 1');
    if (r <= d/alpha0)
        xp = x; yp = plus(y);
        %ev = 1;
    elseif (r <= 2*d/alpha0)
        xp = plus(x); yp = y;
        %ev = 2;
    elseif (r <= 3*d/alpha0)
        xp = x; yp = minus(y);
        %ev = 3;
    elseif (r <= 4*d/alpha0)
        xp = minus(x); yp=y;
        %ev = 4;
    elseif (r <= 5*d/alpha0)
        xp = plus(x); yp = minus(y);  
        %ev = 5;
    else
        xp = minus(x); yp = plus(y);
        %ev = 6;
    end
    % update X
    X(x,y) = X(x,y) - 1; % move particle from (x,y)
    X(xp, yp) = X(xp, yp) + 1;  % to (xp, yp)
elseif (r <= (6*d+k)/alpha0) %production
    %ev = 7;
    X(1,1) = X(1,1) + 1;
else %degradation
    %ev = 8;
    X(x,y) = X(x,y) - 1;
end
%choice(xp,yp) = choice(xp,yp) + 1;
%events(ev) = events(ev) + 1;

% Track progress
    if mod(i, 10^4)==0
        disp(i)
        disp(sum(sum(X)));
    end
%}
end

%% display concentration on lattice 
figure();
imagesc(X);
h = colorbar;
%{ 
%Figures for troubleshooting
%%
figure();
imagesc(choice);
%%
figure();
plot(1:8, events, 'o');
%}
%% Calculate signaling molecule concentration as function of distance to cell
% Count the number of molecules at some given distance from the source
[dist, pos] = init_dist_hex(L, L);
distvals = unique(round(dist(1,:), 6)); 
nvals = zeros(1, numel(distvals)); % number of sites with given distance
cvals = zeros(1, numel(distvals)); % concentration of signaling molecule

for i=1:numel(distvals)
    idx = find(round(dist(1,:), 6) == distvals(i));
    nvals(i) = numel(idx);
    cvals(i) = sum(X(idx));
end

% normalised concentration
cvalsN = cvals./nvals./cvals(1); %normalise by c(0) = 1;

figure();
hold on
%plot(distvals, cvalsN, 'o-');
semilogy(distvals, cvalsN, 'o-');
%ylim([200 350]);
%cR = @(r) exp(-r)/(r+1);
%plot(distvals(2:end), arrayfun(cR,distvals(2:end)));

cvalsN1 = cvalsN; 

%% Store calculated concentration
fname = strrep(sprintf('L%d_k_10e%.1f_gamma_10e%.1f_d_10e%.1f_X0_%d_nsteps%d',...
    L, log(k)/log(10), log(gamma)/log(10), log(d)/log(10), X0, nsteps), '.', 'p');
out_file = fullfile(pwd, 'data', 'hex_lattice_diffusion', fname);
%save(out_file, 'distvals', 'cvalsN', 'cvals', 'nvals', 'L', 'k', 'gamma', 'd', 'nsteps');

%% Load previous data
%{
nsteps = 10^6;
fname = strrep(sprintf('L%d_k_10e%.1f_gamma_10e%.1f_d_10e%.1f_X0_%d_nsteps%d',...
    L, log(k)/log(10), log(gamma)/log(10), log(d)/log(10), X0, nsteps), '.', 'p');
in_file = fullfile(pwd, 'data', 'hex_lattice_diffusion', fname);
load(in_file, 'distvals', 'cvalsN', 'cvals', 'nvals', 'L', 'k', 'gamma', 'd', 'nsteps');

figure();
hold on
plot(distvals, cvalsN);
plot(distvals, cvalsN1);
%plot(distvals, cvalsN- cvalsN1);

%}
%% Calculate matrix of signalling strengths
%{
dist2 = round(dist, 6);
distvals = unique(dist2(1,:));

dist_idx = zeros(L^2,L^2);
M = zeros(L^2, L^2);
for i=1:L^2
    for j=1:L^2
        idx = find(dist2(i,j) == distvals, 1);
        dist_idx(i,j) = idx;
        M(i,j) = cvalsN(idx);
    end
end
%}