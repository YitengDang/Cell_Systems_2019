% This calculates using the integral approximation:
% 1 - The minimum number of ON neighbors required to keep an ON cell ON
% 2 - The maximum number of ON cells possible to keep an OFF cell OFF

clear variables
close all

K = linspace(1, 20, 1000); % K's to be tested
n = 0:6;
% Parameters of the system
gridsize = 11;
N = gridsize^2;
a0 = 0.5;
Rcell = 0.2*a0;
save = 1;
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
p = 0.3; % ratio of ON cells
p = ceil(p*N)/N;
for i = 1:numel(n)
    Son_on(:,i) = (K-1-fN)./(1 + (beta1-beta2/N)*n(i) + (p-1/N)*beta2) + 1;
    Son_off(i,:) = (K-1-fN)./((beta1-beta2/N)*n(i) + p*beta2) + 1;
end

n_var = 6*p*(1-p)*(N-6)/(N-1);

figure(1)
plot(K, Son_on)
ylim([1 30])
hold on
%plot(K, K, 'k--')
plot([fN+1 fN+1], [1 30], 'k--')
plot(K, (K-1)/fN, 'b--')
plot(K, K/(fN+1), 'k--')
plot(K, K-fN, 'r--')
hold off
legend({'0', '1', '2', '3', '4', '5', '6'}, 'Location', 'northwest')
set(gca,'ydir','normal', 'FontSize', 16)
xlabel('K', 'FontSize', 18)
ylabel('S_{ON}', 'FontSize', 18)
title(sprintf('ON cells, p = %.1f', p), 'FontSize', 16)

figure(2)
plot(K, Son_off)
ylim([1 30])
hold on
%plot(K, K, 'k--')
plot([fN+1 fN+1], [1 30], 'k--')
plot(K, (K-1)/fN, 'b--')
plot(K, K/(fN+1), 'k--')
plot(K, K-fN, 'r--')
hold off
legend({'0', '1', '2', '3', '4', '5', '6'}, 'Location', 'southeast')
set(gca,'ydir','normal', 'FontSize', 16)
xlabel('K', 'FontSize', 18)
ylabel('S_{ON}', 'FontSize', 18)
title(sprintf('OFF cells, p = %.1f', p), 'FontSize', 16)