% Runs a single outcome of the Monte Carlo simulation using the Langevin
% approach. It considers dp and dtheta are dependent and estimates the
% covariance of them.

% Update from 01-Nov-2016: Disregarding covariance is working better

clear variables
close all

% Parameters of the system
gridsize = 11; % size of re hexagonal grid
a0 = 1.5;
Rcell = 0.2*a0;
Con = 8;
K = 16;
alpha = 2; % noise

% Initial state (the run is for one state)
p_ini = 0.6;
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
xiN = sum(aux(1,2:end));

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
    p_minus = p*(1-Ponon);
    p_plus = (1-p)*(1-Poffoff);
    
    dp = p_plus - p_minus;
    
    var_p_plus = p_plus*Poffoff/N;
    var_p_minus = p_minus*Ponon/N;
    vardp = var_p_plus + var_p_minus;
    
    % Calculate the average and variance of dtheta. Take care that the
    % function has a problem when P(off->off) and/or P(on->on) are close to
    % unity.
    if 1 - Poffoff < 1e-16
        y_plus = muoff;
        var_y_plus = sigmaoff^2;
    else
        y_plus = muoff + sigmaoff*phi(zoff)/(1 - Poffoff);
        if Poffoff < 1e-16
            var_y_plus = 0;
        else
            var_y_plus = sigmaoff^2*(1 + zoff*phi(zoff)/(1 - Poffoff) ...
                - (phi(zoff)/(1 - Poffoff))^2);
        end
    end
    
    if 1 - Ponon < 1e-16
        y_minus = muon;
        var_y_minus = sigmaon^2;
    else
        y_minus = muon + sigmaon*phi(zon)/(1 - Ponon);
        if Ponon < 1e-16
            var_y_minus = 0;
        else
            var_y_minus = sigmaon^2*(1 - zon*phi(zon)/(1 - Ponon) ...
                - (phi(zon)/(1 - Ponon))^2);
        end
    end
    
    sigma_minus = p_minus^2*var_y_minus + ...
        y_minus^2*var_p_minus + var_y_minus*var_p_minus;
    
    sigma_plus = p_plus^2*var_y_plus + ...
        y_plus^2*var_p_plus + var_y_plus*var_p_plus;
    
    dtheta = 8/(Con-1)*(p_plus*y_plus - p_minus*y_minus) ...
        + 4*fN*(dp^2+vardp - dp*(Con+1)/(Con-1));
    
    vardtheta = 64/(Con-1)^2*(sigma_plus + sigma_minus) + ...
        16*fN^2*vardp*((Con+1)^2/(Con-1)^2 + 4*dp^2 + 2*vardp);
    
    % Estimate the covariance of p and theta
    covdpdtheta = 8/(Con-1)*(y_plus*(var_p_plus + p_plus*dp) ...
        + y_minus*(var_p_minus - p_minus*dp)) ...
        + 4*fN*(dp^3 + 3*dp*vardp - ...
        (Con+1)/(Con-1)*(vardp + dp^2)) - dp*dtheta;
    
    % Evolve the state one time step
    if vardp <=0
        % If vardp is zero then there is no noise on p
        theta = theta + normrnd(dtheta, sqrt(vardtheta));
        p = p + dp;
    else
        cov_d = [vardp covdpdtheta; covdpdtheta vardtheta];
        mean_d = [dp dtheta];
        
        dx = mvnrnd(mean_d, cov_d);
        
        p = p + dx(1);
        % Check if p moved in the step
        if dx(1) < min_dp
            disp('nomove')
        end
        theta = theta + dx(2);
    end
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

