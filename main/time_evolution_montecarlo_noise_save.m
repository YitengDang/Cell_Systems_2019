% Time evolution of a system with noise and without visualization

close all
clear all
warning off

%data_path = 'C:\Users\eduardopavinat\Dropbox\Matlab codes\data_onecelltype_entropy';

% Parameters of the system
gridsize = 15;
N = gridsize^2;
a0 = 0.5;
Rcell = 0.2*a0;
p = 0.8;
iniON = round(p*N);
I_ini = 0;
n_run = 50;
tmax = 100;
% parameters
Con = 8;
K = 15;
alpha = 0.8;

% use hexagonal lattice
[dist, pos] = init_dist_hex(gridsize, gridsize);

dist_vec = a0*dist(1,:);
r = dist_vec(dist_vec>0); % exclude self influence
fN = sum(sinh(Rcell)*sum(exp(Rcell-r)./r)); % calculate signaling strength
gN = sum(sum((sinh(Rcell)*exp(Rcell-r)./r).^2)); % calculate signaling strength

hfunc = @(p, I) -0.5*(Con-1)*(1 + 4*fN.*p.*(1-p).*I + fN*(2*p-1).^2) ...
    -(2*p-1).*(0.5*(Con+1)*(1+fN) - K);

for tests = 1:n_run
    disp(tests)
    I_t = I_ini + normrnd(0,0.01);
    p = iniON/N;
    theta = fN*(4*I_t*p*(1-p)+(2*p-1)^2);
    t = 0;
    I = [];
    Non = [];
    mom = [];
    Non(end+1) = round(p*N);
    I(end+1) = I_t;
    mom(end+1) = hfunc(p,I_t);
    for i = 1:tmax
        [theta, p, ~] = update_montecarlo(theta, p, N, Con, K, fN, gN, alpha);
        t = t+1;
        I_t = calc_I(p, theta, fN);
        Non(end+1) = round(p*N);
        I(end+1) = I_t;
        mom(end+1) = hfunc(p,I_t);
    end
    
    fname_str = sprintf('N%d_n%d_neq_%d_a0%d_K_%d_Son_%d_t_%d_noise_%d_montecarlo', ...
        N, iniON, Non(end), 10*a0, K, Con, t, 10*alpha);
    i = 1;
    fname = fullfile(pwd, 'data', 'dynamics', 'noise_MC',...
        strcat(fname_str,'-v',int2str(i),'.mat'));
    while exist(fname, 'file') == 2
        i=i+1;
        fname = fullfile(pwd, 'data', 'dynamics', 'noise_MC',...
            strcat(fname_str,'-v',int2str(i),'.mat'));
    end
    
    save(fname)
end
