% Simulates Langevin equation based on gradient of Hamiltonian with random
% Gaussian noise
clear all
close all
clc
%% Parameters
gridsize = 11;
N = gridsize^2;
a0 = 0.5;
Rcell = 0.2*a0;
K = 16;
Con = 8;

nruns = 10; % number of simulations
maxsteps = 100; % number of steps per simulation
delta = 0.01; % step size in Langevin equation
noise = 2; % Lagevin noise amplitude
alpha = 0; % sensing noise
% initial conditions
p_ini = 0.5;
I_ini = 0;

% load fN
[dist, pos] = init_dist_hex(gridsize, gridsize);
dist_vec = dist(1,:);
r = a0*dist_vec(dist_vec>0); % exclude self influence
fN = sum(sinh(Rcell)*sum(exp(Rcell-r)./r)); % calculate signaling strength sphere
gN = sum(sum((sinh(Rcell)*exp(Rcell-r)./r).^2)); % calculate signaling strength

%% Simulation
% define gradient
h = @(p, I) -0.5*(Con-1)*(1 + 4*fN.*p.*(1-p).*I + fN*(2*p-1).^2) ...
    -(2*p-1).*(0.5*(Con+1)*(1+fN) - K);
delh_delp = @(p,I) (Con-1)*2*fN*(I-1).*(2*p - 1) -(Con+1)*(fN+1)+2*K;
delh_delI = @(p,I) 2*fN*(Con - 1)*p.*(p - 1);

% trajectory variables
p = zeros(nruns, maxsteps+1);
I = zeros(nruns, maxsteps+1);
p(:,1) = p_ini;
I(:,1) = I_ini;

% run simulation
for i=1:nruns
    for j=1:maxsteps 
        pt = p(i,j); 
        It = I(i,j);
        [peq, p_lim] = revision_update_Peq(pt, It, N, Con, K, fN, gN, alpha); %
        if rand() < peq % system in equilibrium?
            if p_lim == 1 || p_lim == -1 %if system goes out of range
                p(i, j:end) = (p_lim+1)/2;
                I(i, j:end) = 0;
            else
                p(i, j:end) = pt;
                I(i, j:end) = It;
            end
            break
        end
        p(i,j+1) = p(i,j) - delh_delp(pt, It)*delta + normrnd(0,noise)*delta;
        I(i,j+1) = I(i,j) - delh_delI(pt, It)*delta + normrnd(0,noise)*delta;
    end
end

%% Plot result
figure();
hold on

% Plot vector field
pv = (0:5:N)/N;
Iv = -1:0.1:1;
[pm, Im] = meshgrid(pv, Iv);
vf_x = -delh_delp(pm, Im);
vf_y = -delh_delI(pm, Im);
quiver(pm, Im, vf_x, vf_y, 'LineWidth', 0.8, 'AutoScaleFactor', 1.2);

% plot trajectories
plot(p', I', 'g-');
%xlim([0 1]);
ylim([-0.1 1.05]);

% plot start and end points
plot(p(:,1), I(:,1), 'go');
plot(p(:,end), I(:,end), 'gx');
