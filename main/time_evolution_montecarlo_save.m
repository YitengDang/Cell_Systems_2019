% Time evolution of a system without noise and without visualization
close all
clear all
warning off
set(0, 'defaulttextinterpreter', 'latex');

% Parameters of the system
gridsize = 15;
N = gridsize^2;
a0 = 0.5;
Rcell = 0.2*a0;
% circuit parameters
Con = 11;
K = 18;
% initial conditions
p0 = 0.5;
%p0_all = (0:5:N)/N;
iniON = round(p0*N);
I0 = 0;

% simulation parameters
n_run = 50;

% use hexagonal lattice
[dist, pos] = init_dist_hex(gridsize, gridsize);

dist_vec = a0*dist(1,:);
r = dist_vec(dist_vec>0); % exclude self influence
fN = sum(sinh(Rcell)*sum(exp(Rcell-r)./r)); % calculate signaling strength
gN = sum(sum((sinh(Rcell)*exp(Rcell-r)./r).^2)); % calculate signaling strength

hfunc = @(p, I) -0.5*(Con-1)*(1 + 4*fN.*p.*(1-p).*I + fN*(2*p-1).^2) ...
    -(2*p-1).*(0.5*(Con+1)*(1+fN) - K);

for run = 1:n_run %numel(p0_all)
    disp(run)
    %p0 = p0_all(run);
    %iniON = round(p0*N);

    I_t = I0;
    p = iniON/N;
    theta = fN*(4*I_t*p*(1-p)+(2*p-1)^2);
    t = 0;
    I = [];
    I_2 = [];
    Non = [];
    h = [];
    Non(end+1) = round(p*N);
    I(end+1) = I_t;
    I_2(end+1) = I_t;
    h(end+1) = hfunc(p,I_t);
    %[theta, p, pe] = update_montecarlo_2(theta, p, N, Con, K, fN, gN, 0);
    [theta, theta_2, p, pe] = update_montecarlo_YD(theta, p, N, Con, K, fN, gN, 0);
    while rand > pe
        t = t+1;
        Non(end+1) = round(p*N);
        I(end+1) = calc_I(p, theta, fN);
        I_2(end+1) = calc_I(p, theta_2, fN);
        h(end+1) = hfunc(p,I_t);
        %[theta, p, pe] = update_montecarlo_2(theta, p, N, Con, K, fN, gN, 0);
        [theta, theta_2, p, pe] = update_montecarlo_YD(theta, p, N, Con, K, fN, gN, 0);
    end
    
    fname_str = strrep(sprintf('N%d_n%d_I0_%.1f_neq_%d_a0_%.2f_K_%d_Con_%d_t_%d_montecarlo_v3_dp_stoch', ...
        N, iniON, I0, Non(end), a0, K, Con, t), '.', 'p');
    i = 1;
    fname = fullfile(pwd, 'data', 'dynamics', 'no_noise_MC',...
        strcat(fname_str,'-v',int2str(i),'.mat'));
    while exist(fname, 'file') == 2
        i=i+1;
        fname = fullfile(pwd,'data', 'dynamics', 'no_noise_MC',...
            strcat(fname_str,'-v',int2str(i),'.mat'));
    end
    
    save(fname)
end
