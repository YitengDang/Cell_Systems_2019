% Time evolution of a system without noise and without visualization
close all
clear all
warning off

% Parameters of the system
gridsize = 11;
N = gridsize^2;
a0 = 0.5;
Rcell = 0.2*a0;
%Con = 14;
%K = 17;

nrange = 0:N;
Klist = [10 10 14 19 16 18]; % (a0=1.5) [3 6 10 15 17 20]; % (a0=0.5) [10 10 14 19 16 18];
Conlist = [5 21 16 14 8 6]; % (a0=1.5) [24 21 21 20 14 14]; % (a0=0.5) [5 21 16 14 8 6];

% Simulation parameters
n_run = 1; % number of runs per value of K, Con, iniON
%p = 0.5;
%iniON = round(p*N);
I_ini = 0;

% use hexagonal lattice
[dist, pos] = init_dist_hex(gridsize, gridsize);
dist_vec = a0*dist(1,:);
r = dist_vec(dist_vec>0); % exclude self influence
fN = sum(sinh(Rcell)*sum(exp(Rcell-r)./r)); % calculate signaling strength
gN = sum(sum((sinh(Rcell)*exp(Rcell-r)./r).^2)); % calculate signaling strength

for i=1:length(Klist)
    K = Klist(i);
    Con = Conlist(i);
    
    % hamiltonian
    hfunc = @(p, I) -0.5*(Con-1)*(1 + 4*fN.*p.*(1-p).*I + fN*(2*p-1).^2) ...
        -(2*p-1).*(0.5*(Con+1)*(1+fN) - K);

    for iniON = nrange
        fprintf('n=%d \n', iniON);
        for tests = 1:n_run
            disp(tests);
            I_t = I_ini;
            p = iniON/N;
            theta = fN*(4*I_t*p*(1-p)+(2*p-1)^2);
            t = 0;
            I = [];
            Non = [];
            mom = [];
            Non(end+1) = round(p*N);
            I(end+1) = I_t;
            mom(end+1) = hfunc(p,I_t);
            [theta, p, pe] = update_montecarlo_2(theta, p, N, Con, K, fN, gN, 0);
            while rand > pe
                t = t+1;
                I_t = calc_I(p, theta, fN);
                Non(end+1) = round(p*N);
                I(end+1) = I_t;
                mom(end+1) = hfunc(p,I_t);
                [theta, p, pe] = update_montecarlo_2(theta, p, N, Con, K, fN, gN, 0);
            end
            
            fname_str = sprintf('N%d_n%d_I0_%d_neq_%d_a0%d_K_%d_Son_%d_t_%d_MonteCarlo', ...
                N, iniON, I_ini, Non(end), 10*a0, K, Con, t);
            i = 1;
            fname = fullfile(pwd, 'rebuttal', 'h_vector_field', 'data', ...
                strcat(fname_str,'-v',int2str(i),'.mat'));
            while exist(fname, 'file') == 2
                i=i+1;
                fname = fullfile(pwd, 'rebuttal', 'h_vector_field', 'data', ...
                    strcat(fname_str,'-v',int2str(i),'.mat'));
            end

            save(fname)
        end
    end
end