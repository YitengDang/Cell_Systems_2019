function theta_max = max_theta(p, a0, Rcell, N, fN)
% Estimate the maximum theta for a fraction of ON cells.
% We use the continous approximation to do so.

G = 2*pi*exp(Rcell)*sinh(Rcell);
Acell = 0.5*sqrt(3)*a0^2; % approximate area of a cell

% Assume there is full organization, all ON cells are within a certain
% distance
Lon = sqrt((p*N-7)*Acell/pi + (1.5*a0)^2);
Loff = sqrt(((1-p)*N-7)*Acell/pi + (1.5*a0)^2);
L = sqrt((N-7)*Acell/pi + (1.5*a0)^2);

theta_max = fN - 2*G*(p*exp(-Lon) + (1-p)*exp(-Loff) - exp(-L));