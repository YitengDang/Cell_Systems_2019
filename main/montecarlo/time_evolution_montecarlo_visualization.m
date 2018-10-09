%% Plot trajectories (simulated, Langevin) on top of a map of Peq in p, I space
clear all
close all
clc
set(0, 'defaulttextinterpreter', 'latex');
%% Parameters
gridsize = 15;
N = gridsize^2;
a0 = 3;
rcell = 0.2;
Rcell = rcell*a0;
K = 11; % A1: 2-4, A01: 5-12, A0: 13-15
Con = 12;
alpha = 0; % noise

% load fN, gN
[dist, pos] = init_dist_hex(gridsize, gridsize);
dist_vec = dist(1,:);
r = a0*dist_vec(dist_vec>0); % exclude self influence
fN = sum(sinh(Rcell)*sum(exp(Rcell-r)./r)); % calculate signaling strengthN
gN = sum(sum((sinh(Rcell)*exp(Rcell-r)./r).^2)); % calculate signaling strength

% Simulation parameters
p0 = 0.5;
I0 = 0;
n_runs = 10;
tmax = 1000;
noSpatialOrder = 0;

% filename for saving
fname_str = strrep(sprintf(...
    'Montecarlo_repr_N%d_a0_%.2f_K_%d_Con_%d_noise_%.1f_p_ini_%.1f_I_ini_%.1f_withoutI_%d_t_%d',...
	N, a0, K, Con, alpha, p0, I0, noSpatialOrder, tmax), '.', 'p');
%% Single simulation 
p_t = [];
%I_t = [];
Theta_t = [];

p = p0;
theta = fN*((2*p-1)^2 + 4*p*(1-p)*I0);

t = 0;
p_t(end+1) = p;
Theta_t(end+1) = theta;
%[theta, p, pe] = update_montecarlo(theta, p, N, Con, K, fN, gN, alpha);
[theta, p, pe] = update_montecarlo_repression(theta, p, N, Con, K, fN, gN, alpha, noSpatialOrder);
while rand > pe && t < tmax
    t = t+1;
    p_t(end+1) = p;
    Theta_t(end+1) = theta;
    %[theta, p, pe] = update_montecarlo(theta, p, N, Con, K, fN, gN, alpha);
    [theta, p, pe] = update_montecarlo_repression(theta, p, N, Con, K, fN, gN, alpha, noSpatialOrder);
end

% save trajectory
qsave = 0;
if qsave
    folder = 'H:\My Documents\Multicellular automaton\figures\one_signal_repression';
    i = 1;
    fname = fullfile(folder, strcat(fname_str,'-v',int2str(i), '.mat'));
    while exist(fname, 'file') == 2
        i=i+1;
        fname = fullfile(folder, ...
            strcat(fname_str,'-v',int2str(i),'.mat'));
    end
    save(fname);
end
%% Plot p, Theta vs time
tmin = 0;
tmax = min(100, t);

h = figure;
hold on
legend_lbl = {'$$p$$'};
plot(tmin:tmax, p_t(tmin+1:tmax+1), 'b-');
if ~noSpatialOrder
    plot(tmin:tmax, Theta_t(tmin+1:tmax+1)/fN, 'r-');
    legend_lbl{end+1} = '$$\Theta/f_N$$';
end
%legend(legend_lbl, 'Interpreter', 'latex');
ylim([0 1]);
set(gca, 'FontSize', 20);
xlabel('t');
ylabel('p');

% save figure
qsave = 0;
if qsave
    folder = 'H:\My Documents\Multicellular automaton\figures\one_signal_repression';
    i = 1;
    fname = fullfile(folder, strcat(fname_str,'-v',int2str(i)));
    while exist(fname, 'file') == 2
        i=i+1;
        fname = fullfile(folder, ...
            strcat(fname_str,'-v',int2str(i)));
    end
    save_figure(h, 10, 8, fname, '.pdf');
end
%% Troubleshoot EOM (update_montecarlo*)
%{
p = p0;
I = 0;

if p <= 0 || p >= 1
    I = 0;
else
    I = (theta - (2*p-1)^2*fN)/4/p/(1-p)/fN;
end
muon = fN*(Con*p + 1 - p + (Con-1)*(1-p)*I);
muoff = fN*(Con*p + 1 - p - (Con-1)*p*I);
kappa = sqrt((Con-1)^2*(gN)*p*(1-p));
sigmaon = real(sqrt(kappa^2+alpha^2));
sigmaoff = real(sqrt(kappa^2+alpha^2));

zoff = (K - 1 - muoff)/sigmaoff;
zon = (K - Con - muon)/sigmaon;
Poffoff = 1-normcdf(zoff);
Ponon = normcdf(zon);

pe = (Ponon^p*Poffoff^(1-p))^N;

% Calculate the average fraction of cells that change state and the
% variance
%     p_minus_mean = p*(1-Ponon);
%     p_plus_mean = (1-p)*(1-Poffoff);
%     
%     var_p_plus = p_plus_mean*Poffoff/N;
%     var_p_minus = p_minus_mean*Ponon/N;
%     
%     if var_p_plus > 0
%         p_plus = normrnd(p_plus_mean, sqrt(var_p_plus));
%     else
%         p_plus = p_plus_mean;
%     end
%     
%     if var_p_minus > 0
%         p_minus = normrnd(p_minus_mean, sqrt(var_p_minus));
%     else
%         p_minus = p_minus_mean;
%     end

p_minus = binornd(round(N*p), 1-Ponon)/N;
p_plus = binornd(round(N*(1-p)), 1-Poffoff)/N;

% Calculate the average and variance of dtheta. Take care that the
% function has a problem when P(off->off) and/or P(on->on) are close to
% unity.
naux = round(p_minus*N);
if naux > 0
    if sigmaon > 0
        zon = (K - Con - muon)/kappa*ones(naux,1);
        Y_minus = mean(muon + ...
            kappa*trandn( zon, (fN*Con-muon)/kappa*ones(naux,1) ) );
    elseif sigmaon==0 % then all cells are ON
        Y_minus = fN*Con;
    end
else
    Y_minus = 0;
end

naux = round(p_plus*N);
if naux > 0
    if sigmaoff > 0
        zoff = (K - 1 - muoff)/kappa*ones(naux,1);
        Y_plus = mean(muoff + ...
            kappa*trandn( (fN-muoff)/kappa*ones(naux,1), zoff ) );
    elseif sigmaoff == 0 % then all cells are OFF
        Y_plus = fN;
    end
else
    Y_plus = 0;
end

% Evolve the state one time step
dp = p_plus - p_minus;
% dp = -(1-Ponon)*p+(1-Poffoff)*(1-p); % alternative (not good)
dtheta1 = 8/(Con-1)*(p_plus*Y_plus - p_minus*Y_minus);
dtheta2 = -4*(Con+1)/(Con-1)*fN*dp;
dtheta3 = 4*fN*dp^2;
dtheta = dtheta1 + dtheta2 + dtheta3;

p_new = p + dp;
p_new = max(0,min(p_new,1));
%{
if theta_new/fN == 1
    p_new = round(p);
else
    p_new = p + dp;
    p_new = round(N*p_new)/N; % To make it a possible value
    p_new = max(0,min(p_new,1));
end
%}

theta_new = theta + dtheta;
theta_new = min(theta_new, fN);
%}
%% Plot on top of Peq map 
%{
% Calculate P_eq
p = (0:N)/N;
I = linspace(-0.3, 1, 250);
[pm, Im] = meshgrid(p,I);

muon = fN*(Con.*pm + 1 - pm + (Con-1).*(1-pm).*Im);
muoff = fN*(Con.*pm + 1 - pm - (Con-1).*pm.*Im);
kappa = sqrt((Con-1)^2*(gN)*pm.*(1-pm));
sigmaon = kappa;
sigmaoff = kappa;

zon = (K - Con - muon)./sigmaon;
zoff = (K - 1 - muoff)./sigmaoff;

Poffoff = normcdf(zoff);
Ponon = 1-normcdf(zon);

pe = (Ponon.^pm .* Poffoff.^(1-pm)).^N;

%% Plot result
h1=figure(idx);
imagesc(p, I, pe);
set(gca, 'YDir', 'Normal');
c = colorbar;
xlabel('p');
ylabel('I');
ylabel(c, 'P_{eq}');
caxis([0 1]);
set(gca, 'FontSize', 24);
set(gcf, 'Units','Inches', 'Position', [0 0 10 8]);

%% Plot gradient vector field
% Multiplicative noise in drift of I
%
%{
% old approach (only N)
Ns = 6^2; %optimal N
fmax = 5; %max amplification
a = Ns*log(fmax);
f = fmax.*exp(-a./N);
f = 1;
%}

% Load fitted surface phi(N, a0)
%fname = 'H:\My Documents\Multicellular automaton\rebuttal\I_correction\I_correction_phi_N_a0_K6_Con15.mat';
%load(fname, 'f');
%phi = f(gridsize, a0);
% phi = maxI Langevin / maxI simulation
% In update rule below, enter 1/phi instead of phi!

% define gradient
delh_delp = @(p,I) (Con-1)*2*fN*(I-1).*(2*p - 1) -(Con+1)*(fN+1)+2*K;
delh_delI = @(p,I) 2*fN*(Con - 1)*p.*(p - 1);

figure(h1);
hold on
pv = (0:gridsize:N)/N;
Iv = -0.3:0.08:1;
[pm2, Im2] = meshgrid(pv, Iv);
vf_x = -delh_delp(pm2, Im2);
vf_y = -delh_delI(pm2, Im2);
quiver(pm2, Im2, vf_x, vf_y, 'LineWidth', 1, 'AutoScaleFactor', 1.2,...
    'Color', [0.5 0.5 0.5]);
%}
