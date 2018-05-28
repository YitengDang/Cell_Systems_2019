clear variables
close all

% parameter of the run
gridsize = 30;
a0 = 0.5;
Rcell = 0.2*a0;
Con = 5;
K = 10;
alpha = 0; % noise

% Initial state
p_ini = 0.6; % fraction of ON cells
Ii = 0; % spatial order parameter

% Place cells in hexagonal grid
N = gridsize^2;
[dist, pos] = init_dist_hex(gridsize, gridsize);

dist_vec = dist(1,:);
r = a0*dist_vec(dist_vec>0); % exclude self influence

% Get the signalling length
fN = sum(sinh(Rcell)*sum(exp(Rcell-r)./r)); % calculate signaling strength sphere
gN = sum((sinh(Rcell)^2)*sum(exp(2*(Rcell-r))./r.^2)); % calculate signaling strength sphere

% Calculate the approximate dynamics of p and H

% calculate hamiltonian given p and I
hfunc = @(p, I) -0.5*(Con-1)*(1 + 4*fN.*p.*(1-p).*I + fN*(2*p-1).^2) ...
    -(2*p-1).*(0.5*(Con+1)*(1+fN) - K);
% auxiliary function
phi = @(x) exp(-0.5*x^2)/sqrt(2*pi);

% initial hamiltonian
hi = hfunc(p_ini, Ii);

% number of steps to run
steps = 100;
% initialize variables
h = zeros(1,steps);
p = zeros(1,steps);
I = zeros(1,steps);
theta = zeros(1,steps);
pe = zeros(1,steps-1);
p(1) = p_ini;
I(1) = Ii;
h(1) = hi;
theta(1) = 4*p_ini*(1-p_ini)*fN*Ii+(2*p_ini-1)^2*fN;
% run for all the steps
for i = 2:steps
    % std of the neighbor's concentration
    sigmap = sqrt((Con-1)^2*gN*p(i-1)*(1-p(i-1))+alpha^2);
    % Calculate P(ON->ON) and P(OFF->OFF)
    zon = Con - K + 0.5*(Con+1)*fN + 0.5*(Con-1)*((2*p(i-1)-1)*fN + 2*fN*(1-p(i-1))*I(i-1));
    zoff = 1 - K + 0.5*(Con+1)*fN + 0.5*(Con-1)*((2*p(i-1)-1)*fN - 2*fN*p(i-1)*I(i-1));
    Poffoff = normcdf(-zoff/sigmap);
    Ponon = normcdf(zon/sigmap);
    % calculate the average change in the fraction of ON cells
    ponoff = p(i-1)*(1-Ponon);
    poffon = (1-p(i-1))*(1-Poffoff);
    dp = poffon-ponoff;
    % calculate the average change in theta
    dtheta = 8/(Con-1)*((K-0.5*(Con+1)*fN)*dp + Con*ponoff - poffon) + ...
        4*fN*(ponoff^2 + poffon^2 - 2*ponoff*poffon);
    % evolve the system
    theta(i) = theta(i-1) + dtheta;
    p(i) = p(i-1) + dp;
    pe(i-1) = Ponon^p(i-1) * Poffoff^(1-p(i-1));
    % Calculate the new spatial order parameter
    eps = 1e-6;
    if p(i-1) < eps || p(i-1) > 1-eps
        I(i) = 0;
    else
        dI = (0.25*dtheta/fN - (2*p(i-1)-1)*(1-I(i-1))*dp)/p(i-1)/(1-p(i-1));
        I(i) = I(i-1) + dI;
    end
    % calculate the new Hamiltonian
    h(i) = hfunc(p(i), I(i));
end

% Plot the fraction of ON cells in time
figure(1)
plot(p)
xlabel('Time (Steps)')
ylabel('Fraction p')

% Plot the spatial order parameter in time
figure(2)
plot(I)
xlabel('Time (Steps)')
ylabel('Spatial order')

% Plot the equilibrium probability in time
figure(3)
plot(pe)
xlabel('Time (Steps)')
ylabel('Equilibrium prob.(p_e^{1/N})')

% Plot the hamiltonian in time
figure(4)
plot(h)
xlabel('Time (Steps)')
ylabel('Hamiltonian (H/N)')

% Plot the order theta in time
figure(5)
plot(theta/fN, '-x')
xlabel('Time (Steps)')
ylabel('\Theta')

% plot the path in p vs. I
figure(6)
plot(p,I, '-x')

