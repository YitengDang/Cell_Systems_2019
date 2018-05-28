% Runs a single outcome of the Monte Carlo simulation using the Langevin
% approach. It considers dp and dtheta are dependent and estimates the
% covariance of them.

% Update from 01-Nov-2016: Disregarding covariance is working better

clear variables
close all

% Parameters of the system
gridsize = 15; % size of re hexagonal grid
a0 = 0.5;
Rcell = 0.2*a0;
Con = 8;
K = 15;
alpha = 0; % noise

% Initial state (the run is for one state)
p_ini = 0.7;
Ii = 0;
N = gridsize^2;

% Minimum change in dp (it is a discrete variable)
min_dp = 1/N;

% initialize the hexagonal grid
[dist, pos] = init_dist_hex(gridsize, gridsize);
dist_vec = dist(1,:);
r = a0*dist_vec(dist_vec>0); % exclude self influence

% Get the signalling length
fN = sum(sinh(Rcell)*sum(exp(Rcell-r)./r)); % calculate signaling strength sphere

% Matrix to calculate the influence of neighbors
M = sinh(Rcell)*exp(Rcell-a0*dist)./(a0*dist);
M(1:N+1:N^2) = 0;
% variables used to calculate the std of neighbors concentration
aux = M*M';
gN = aux(1,1);

% Calculate the approximate dynamics of p and H

% Function calculate the Hamiltonian h given the state
hfunc = @(p, I) -0.5*(Con-1)*(1 + 4*fN.*p.*(1-p).*I + fN*(2*p-1).^2) ...
    -(2*p-1).*(0.5*(Con+1)*(1+fN) - K);
% Auxiliary function
phi = @(x) exp(-0.5*x^2)/sqrt(2*pi);

% Initial hamiltonian
hi = hfunc(p_ini, Ii);

% Initialize the variables
pv(1) = p_ini;
Iv(1) = Ii;
hv(1) = hi;
thetav(1) = 4*p_ini*(1-p_ini)*fN*Ii+(2*p_ini-1)^2*fN;
cont = true; % this controls if the loop should continue
i = 1;
while cont
    i = i+1;
    % Get the last state
    p = pv(i-1);
    I = Iv(i-1);
    theta = thetav(i-1);
    % Calculate the average concentration due to neighbors
    muon = fN*(Con*p + 1 - p + (Con-1)*(1-p)*I);
    muoff = fN*(Con*p + 1 - p - (Con-1)*p*I);
    % Calculate the standard deviation
% old version:     kappa = (Con-1)^2*(gN)*p*(1-p)/ ...
%         (p*muon + (1-p)*muoff);
% old version:     sigmaon = sqrt(max(kappa*muon,0)+alpha^2);
% old version:    sigmaoff = sqrt(max(kappa*muoff,0)+alpha^2);
    sigmaon = sqrt((Con-1)^2*(gN)*p*(1-p)+alpha^2);
    sigmaoff = sqrt((Con-1)^2*(gN)*p*(1-p)+alpha^2);
    
    % Calculate P(off->off) and P(on->on)
    zon = (K - Con - muon)/sigmaon;
    zoff = (K - 1 - muoff)/sigmaoff;
    Poffoff = normcdf(zoff);
    Ponon = 1-normcdf(zon);
    % Calculate the equilibrium probability and test if it should continue
    pe(i-1) = (Ponon^p * Poffoff^(1-p))^N;
    if rand < pe(i-1)
       cont = false;
       continue 
    end
    % Calculate the average fraction of cells that change state and the
    % variance
    p_minus_mean = p*(1-Ponon);
    p_plus_mean = (1-p)*(1-Poffoff);
    
    var_p_plus = p_plus_mean*Poffoff/N;
    var_p_minus = p_minus_mean*Ponon/N;
    
    if var_p_plus > 0
        p_plus = normrnd(p_plus_mean, sqrt(var_p_plus));
    else
        p_plus = p_plus_mean;
    end
    
    if var_p_minus > 0
        p_minus = normrnd(p_minus_mean, sqrt(var_p_minus));
    else
        p_minus = p_minus_mean;
    end    
    
    % Calculate the average and variance of dtheta. Take care that the
    % function has a problem when P(off->off) and/or P(on->on) are close to
    % unity.
    if isreal(sigmaon) && sigmaon > 0
        Y_minus = muon + sigmaon*trandn(-inf, zon);
    else
        Y_minus = 0;
    end
    
    if isreal(sigmaoff) && sigmaoff > 0
        Y_plus = muoff + sigmaoff*trandn(zoff, inf);
    else
        Y_plus = 0;
    end

    % Evolve the state one time step
    dp = p_plus - p_minus;
    dtheta = 8/(Con-1)*(p_plus*Y_plus - p_minus*Y_minus) - ...
        4*(Con+1)/(Con-1)*fN*(p_plus - p_minus) + 4*fN*dp^2;

    p = p + dp;
    theta = theta + dtheta;

    % adjust p and theta within their boundaries
    pv(i) = max(0,min(p,1));
    thetav(i) = min(theta, fN);
    
    % Calculate I, taking care of the fact that it is not analytical at p=0
    eps = 1e-6;
    if pv(i) < eps || pv(i) > 1-eps
        Iv(i) = 0;
    else
        Iv(i) = (thetav(i) - (2*pv(i)-1)^2*fN)/4/pv(i)/(1-pv(i))/fN;
        %fprintf('I: %.4f, p: %.4f, theta: %.4f \n',I,p,theta);
        %I(i) = I(i-1) + dI;
    end
    % Calculate the hamiltonian of the new state
    hv(i) = hfunc(pv(i), Iv(i));
end

% Plot the evolution of the fraction of ON cells
figure(1)
plot(pv)
xlabel('Time (Steps)')
ylabel('Fraction p')

% Plot the evolution of the spatial order
figure(2)
plot(Iv)
xlabel('Time (Steps)')
ylabel('Spatial order')

% Plot the evolution of the equilibrium probability
figure(3)
plot(pe)
xlabel('Time (Steps)')
ylabel('Equilibrium prob.(p_e)')

% Plot the evolution of the Hamiltonian
figure(4)
plot(hv)
xlabel('Time (Steps)')
ylabel('Hamiltonian (H/N)')

% Plot the evolution of the order parameter theta
figure(5)
plot(thetav/fN, '-x')
xlabel('Time (Steps)')
ylabel('\Theta')

% Plot the evolution of the spatial order parameter
figure(6)
plot(pv,Iv, '-x')

