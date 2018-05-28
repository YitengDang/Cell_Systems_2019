function [fN, Imax] = max_moranI(a0, Rcell, N, p)
% Estimate the maximum spatial order parameter for a fraction of ON cells.
% We use the continous approximation to do so.

Acell = 0.5*sqrt(3)*a0^2; % approximate area of a cell
G = 2*pi*exp(Rcell)*sinh(Rcell)/Acell;

L0 = 1.5*a0; % Distance from where to assume the continuous approach

% Assume there is full organization, all ON cells are within a certain
% distance
Lon = sqrt(Acell*(p*N-7)/pi + L0^2);
Loff = sqrt(Acell*((1-p)*N-7)/pi + L0^2);
L = sqrt(Acell*(N-7)/pi + L0^2);

fa0 = exp(Rcell-a0)*sinh(Rcell)/a0;
fN = 6*fa0 + G*(exp(-L0) - exp(-L));

% Estimate the maximum I
Imax = (fN-2*G*(p.*exp(-Lon)+(1-p).*exp(-Loff)-exp(-L))-fN*(2*p-1).^2)./ ...
    (4*fN*p.*(1-p));
Imax(or(p==0,p==1)) = 0;