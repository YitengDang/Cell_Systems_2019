% Plot two maps with the number of ON neighbors required to keep cells ON or
% OFF using the nearest neighbor approximation

clear variables
close all

% Parameters of the system
gridsize = 11;
N = gridsize^2;
a0 = 1;
R = 0.2*a0;

Acell = 0.5*sqrt(3)*a0^2;
alpha = 2*pi*R*exp(R)/Acell;
L = sqrt(Acell*N/pi);

Son_vec = linspace(1,30,300);
K_vec = linspace(1,20,200);

[K, Son] = meshgrid(K_vec, Son_vec);

pON = ((K-Son)/alpha/(exp(-a0)-exp(-L)) - 1)./(Son-1) + 1/N;

pON(pON<0 | pON>1) = -1;

figure(1)
imagesc(K_vec, Son_vec, pON)
colorbar;
hold on
plot(K_vec, K_vec, 'r--')
hold off
set(gca, 'ydir', 'normal')

pOFF = ((K-1)/alpha/(exp(-a0)-exp(-L)) - 1)./(Son-1);

pOFF(pOFF<0 | pOFF>1) = -1;

figure(2)
imagesc(K_vec, Son_vec, pOFF)
colorbar;
hold on
plot(K_vec, K_vec, 'r--')
hold off
set(gca, 'ydir', 'normal')
