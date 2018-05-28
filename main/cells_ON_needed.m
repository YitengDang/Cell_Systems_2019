% Plot two maps with the fractionof ON cells required to keep cells ON or
% OFF using the nearest neighbor approximation. This assumes that the
% nearest neighbor approximation is valid and the remaining of the
% population can be approximated by a mean field approach

clear variables
close all

% Parameters of the system
gridsize = 11;
N = gridsize^2;
a0 = 1;
R = 0.2*a0;

Acell = 0.5*sqrt(3)*a0^2; % Area of a cell
alpha = 2*pi*R*exp(R)/Acell; % Constant of the system
L = sqrt(Acell*N/pi); % Maximum size of the system

% Son and K where the calculation is performed
Son_vec = linspace(1,30,300);
K_vec = linspace(1,20,200);

[K, Son] = meshgrid(K_vec, Son_vec);

% Estimate the number of ON cells needed to turn OFF
pON = ((K-Son)/alpha/(exp(-a0)-exp(-L)) - 1)./(Son-1) + 1/N;

% -1 indicates no valid solution
pON(pON<0 | pON>1) = -1;

% plot the result
figure(1)
imagesc(K_vec, Son_vec, pON)
colorbar;
hold on
plot(K_vec, K_vec, 'r--')
hold off
set(gca, 'ydir', 'normal')

% Estimate the number of OFF cells needed to stay OFF
pOFF = ((K-1)/alpha/(exp(-a0)-exp(-L)) - 1)./(Son-1);

% -1 indicates no valid solution
pOFF(pOFF<0 | pOFF>1) = -1;

% plot the result
figure(2)
imagesc(K_vec, Son_vec, pOFF)
colorbar;
hold on
plot(K_vec, K_vec, 'r--')
hold off
set(gca, 'ydir', 'normal')
