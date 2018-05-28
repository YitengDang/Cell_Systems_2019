% Numerically computes the master equation based on the update rule of the
% system
close all
clear all
warning off

% Lattice parameters
gridsize = 11;
N = gridsize^2;
a0 = 5;
Rcell = 0.2*a0;

% circuit parameters
K = 8;
Con = 16;
hill = 2; % Hill coefficient
prec = 5; %precision

% Initial configuration
initialID = 'uniform'; %'uniform' 'allON' 'singleOFF' 'fixedp'

% simulation parameters
tmax = 50;

% calculate fN
[dist, pos] = init_dist_hex(gridsize, gridsize);
dist = round(dist, 5);
dist_vec = a0*dist(1,:);
r = dist_vec(dist_vec>0); % exclude self influence
fN = sum(sinh(Rcell)*sum(exp(Rcell-r)./r)); % calculate signaling strength
%%
ntrials = 10^4;
nbins = 100;
prob = zeros(nbins, tmax+1); % probability distribution p(x,t)
prob(:, 1) = 1/nbins;

% update probability
for t=1:tmax
    for trial = 1:ntrials
        r = rand(N,1);
        x = zeros(N,1); % values of X_i
        cum_prob = cumsum(prob(:, t));
        for i=1:N
            idx = find(r(i) < cum_prob, 1);
            x(i) = idx/nbins;
        end
        y = (1+fN)*sum((Con-1)*x + 1)/N; %sensed concentration
        x_new = y^hill/(K^hill + y^hill); %new state
        bin_idx = floor(x_new*nbins);
        prob(bin_idx, t+1) = prob(bin_idx, t+1) + 1/ntrials;
    end
end

pdfunc = prob*nbins; % probability density from discrete pdf
%% Plot pdf at given time
t = 2;
figure();
binc = (0:1/nbins:1-1/nbins) + 1/nbins/2; %bincenters
plot(binc, pdfunc(:, t+1), 'LineWidth', 2);
xlabel('x');
ylabel(sprintf('P(x,t=%d)', t));
set(gca,'FontSize', 24);
set(gcf, 'Units', 'Inches', 'Position', [0.5 0.5 10 8]);
xlim([0 1]);
ylim([0 max(pdfunc(:, t+1))+1]);

%% Evolution of mean and variance of pdf
binc = (0:1/nbins:1-1/nbins) + 1/nbins/2; %bincenters
xvar = zeros(tmax+1, 1);
xstd = zeros(tmax+1, 1);
for t=0:tmax
    xmean(t+1) = sum( binc.*prob(:, t+1)' );
    xvar(t+1) = sum( (binc.^2).*prob(:, t+1)' ) - xmean(t+1)^2;
end

figure();
hold on
plot(0:tmax, xmean, 'o-');
plot(0:tmax, xvar, 'o-');
ylim([0 1]);

%%
MExmean = xmean;
MExvar = xvar;
MEprob = prob;
clearvars -except MExmean MExvar MEprob binc