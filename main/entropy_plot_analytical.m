close all
clear all
warning off

% Calculate the analytical entropy map for a set of parameters and plot the
% data. It doesn't calculate the Monte Carlo simulation for the entropy.

% Parameters of the system
gridsize = 11;
N = gridsize^2;
a0 = 1;
Rcell = 0.2*a0;

% use hexagonal lattice
[dist, pos] = init_dist_hex(gridsize, gridsize);
dist_vec = dist(1,:);

% --- Entropy as function of Son (fixed K) --- %
K = 10;
Son = linspace(1,40,100)';

% Calculate analytical entropy
omega = zeros(size(Son));
for i=1:size(Son)
    omega(i) = entropy_eq_sphere_2(dist_vec, Son(i), K, a0, Rcell);
end

% Plot the result
figure(4)
plot(Son,omega, 'r-')
xlabel('S_{ON}')
ylabel('Entropy')
%%
% --- Entropy as function of K (Son fixed) --- %
Son = 10;
K = linspace(1,20,100)';

% Calculate analytical entropy
omega = zeros(size(K));
kmax = zeros(size(K));
for i=1:size(K)
    [omega(i),omegak(i,:)] = entropy_eq_sphere_2(dist_vec, Son, K(i), a0, Rcell);
    kmax(i) = find(omegak(i,:) == max(omegak(i,:)),1)-1;
end

% Plot the result
figure(1)
plot(K,log(2^N*omega))
xlabel('K')
ylabel('Entropy')
%%

% Plot Kmax
figure(3)
plotyy(K,kmax, K, omega)
ylabel('k_{max}')

% Plot the result
figure(2)
hold on
plot(0:N,omegak(1,:))
hold off
xlabel('k')
ylabel('\Omega_k')

