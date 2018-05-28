% Estimate the entropy map based on the nearest neighbor approach
% This truns out to be similar to the mean field approach, which we found
% exact
clear variables
close all
warning off

% Parameters of the system
gridsize = 11;
N = gridsize^2;
a0 = 1;
Rcell = 0.2*a0;
K_vec = linspace(1,20,100);
Son_vec = linspace(1,30,100);
save = 0;
fig_file = 'sphere_a0_2_regions';

% hexagonal lattice

L0 = 1.5*a0;
N0 = 7;

Acell = 0.5*sqrt(3)*a0^2;
alpha = 2*pi*Rcell*exp(Rcell)/Acell;
L = sqrt(Acell*(N-N0)/pi + L0^2);

beta1 = Rcell*exp(Rcell-a0)/a0;
beta2 = alpha*(exp(-L0)-exp(-L));
fN = 6*beta1 + beta2;

[K, Son] = meshgrid(K_vec, Son_vec);
omega = zeros(size(K));

for Non = 0:N
    probs = hygepdf(0:6,121,60,6);
    p = Non/N;
    numON = (((K-Son-fN)./(Son - 1) - (p-1/N)*beta2)./(beta1-beta2/N));
    numOFF = (((K-1-fN)./(Son - 1) - p*beta2)./(beta1-beta2/N));
    
    probON = zeros(size(numON));
    probOFF = zeros(size(numOFF));
    for i = 0:min(6, Non)
        % Calculate probability P(Y>K)
        idxON = numON <= i;
        probON(idxON) = probON(idxON) + probs(i+1);
        % Calculate probability P(Y<K)
        idxOFF = numOFF >= i;
        probOFF(idxOFF) = probOFF(idxOFF) + probs(i+1);
    end

    omegak = (probON.^Non).*(probOFF.^(N-Non));
    omega = omega + omegak*nchoosek(N, Non);
end

omega(omega < 1) = 1;
figure(1)
imagesc(K_vec, Son_vec, log(omega))
set(gca, 'ydir', 'normal')
colorbar

figure(2)
plot(Son_vec, log(omega(K<15.1 & K>15)))
        