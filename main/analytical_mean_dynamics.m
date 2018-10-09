% Runs a single outcome of the Monte Carlo simulation using the Langevin
% approach. It considers dp and dtheta as independent variables.

clear variables
close all

% Parameters of the system
gridsize = 25; % size of re hexagonal grid
a0 = 1.5;
Rcell = 0.2*a0;
Con = 35;
K = 35;
alpha = 0; % noise
% Initial state (the run is for one state)
p_ini = 0.5;
Ii = 0.5;
% Number of steps to run
steps = 50;

% initialize the hexagonal grid
N = gridsize^2;
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
xiN = sum(aux(1,2:end));

% Calculate the approximate dynamics of p and H

% Function calculate the Hamiltonian h given the state
hfunc = @(p, I) -0.5*(Con-1)*(1 + 4*fN.*p.*(1-p).*I + fN*(2*p-1).^2) ...
    -(2*p-1).*(0.5*(Con+1)*(1+fN) - K);

% Auxiliary function
phi = @(x) exp(-0.5*x^2)/sqrt(2*pi);

% Initial haniltonian
hi = hfunc(p_ini, Ii);

% Initialize the variables
h = zeros(1,steps);
p = zeros(1,steps);
theta = zeros(1, steps);
I = zeros(1,steps);
pe = zeros(1,steps-1);
p(1) = p_ini;
I(1) = Ii;
h(1) = hi;
theta(1) = 4*p_ini*(1-p_ini)*fN*Ii+(2*p_ini-1)^2*fN;
% test registers the std of dp at each step
test = zeros(1, steps-1);
for i = 2:steps
    % Calculate the average concentration due to neighbors
    muon = fN*(Con*p(i-1) + 1 - p(i-1) + (Con-1)*(1-p(i-1))*I(i-1));
    muoff = fN*(Con*p(i-1) + 1 - p(i-1) - (Con-1)*p(i-1)*I(i-1));
    % Calculate the standard deviation
    kappa = (Con-1)^2*(gN)*p(i-1)*(1-p(i-1))/ ...
        (p(i-1)*muon + (1-p(i-1))*muoff);
% old version: kappa = (Con-1)^2*(gN+xiN*I(i-1)-fN^2*I(i-1)^2)*p(i-1)*(1-p(i-1));
    kappa = max(kappa,0);
    sigmaon = sqrt(kappa*muon+alpha^2);
    sigmaoff = sqrt(kappa*muoff+alpha^2);
% old version: sigmaon = sqrt(kappa + alpha^2);
% old version: sigmaoff = sigmaon;
    
    % Calculate P(off->off) and P(on->on)
    zoff = (K - 1 - muoff)/sigmaoff;
    zon = (K - Con - muon)/sigmaon;
    Poffoff = normcdf(zoff);
    Ponon = 1-normcdf(zon);
    % Calculate the average fraction of cells that change state and the
    % variance
    ponoff = p(i-1)*(1-Ponon);
    poffon = (1-p(i-1))*(1-Poffoff);
    dp = poffon-ponoff;
    vardp = (poffon*Poffoff + ponoff*Ponon)/N;
    test(i-1) = sqrt(vardp);
    
    % Calculate the average change in the order theta
    dtheta = 8/(Con-1)*(poffon*muoff - ponoff*muon + ...
        sigmaoff*(1-p(i-1))*phi(zoff)+ sigmaon*p(i-1)*phi(zon)) ...
        - 4*fN*(Con+1)/(Con-1)*dp + 4*fN*dp^2;
    
    % Calculate the variance of theta
    vartheta = 64/(Con-1)^2*((1-p(i-1))^2*sigmaoff^2* ...
        ((1-Poffoff)^2+zoff*phi(zoff)*(1-Poffoff) - (phi(zoff))^2) ...
        + (p(i-1))^2*sigmaon^2* ...
        ((1-Ponon)^2-zon*phi(zon)*(1-Ponon) - (phi(zon))^2)) + 16*fN^2*dp^4;
    
    % Evolve p and theta
    p(i) = p(i-1) + dp;
    theta(i) = theta(i-1) + dtheta;
    theta(i) = min(theta(i), fN);
    % Calculate the probability of equilibrium of the old state
    pe(i-1) = (Ponon^p(i-1) * Poffoff^(1-p(i-1)))^N;
    % Calculate I, taking care of the fact that it is not analytical at p=0
    eps = 1e-6;
    if p(i) < eps || p(i) > 1-eps
        I(i) = 0;
    else
        I(i) = (theta(i) - (2*p(i)-1)^2*fN)/4/p(i)/(1-p(i))/fN;
        fprintf('I: %.4f, p: %.4f, theta: %.4f \n',I(i),p(i),theta(i));
        %I(i) = I(i-1) + dI;
    end
    % This is to check if there is problem in the calculation of I
    if isnan(I(i))
        disp(kappa)
    end
    % Calculate the hamiltonian of the new state
    h(i) = hfunc(p(i), I(i));
end
%{
% Plot the evolution of the fraction of ON cells
figure(1)
plot(p)
xlabel('Time (Steps)')
ylabel('Fraction p')

% Plot the evolution of the spatial order
figure(2)
plot(I)
xlabel('Time (Steps)')
ylabel('Spatial order')
ylim([-0.05 1])

% Plot the evolution of the equilibrium probability
figure(3)
plot(pe)
xlabel('Time (Steps)')
ylabel('Equilibrium prob.(p_e)')

% Plot the evolution of the Hamiltonian
figure(4)
plot(h)
xlabel('Time (Steps)')
ylabel('Hamiltonian (H/N)')
%}
% Plot the evolution of the order parameter theta
figure(5)
plot(theta/fN, '-x')
xlabel('Time (Steps)')
ylabel('\Theta')
ylim([-0.2 1]);
%{
% Plot the evolution of the spatial order parameter
figure(6)
plot(p,I, '-x')
ylim([-0.05 0.8])

% Plot the variance of dp
figure(7)
plot(test)
%}
