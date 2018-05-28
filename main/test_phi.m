close all
clear all
warning off

% Test for the calculation of the equilibrium probability

% Parameters of the system
gridsize = 11;
N = gridsize^2;
a0 = 0.5;
Rcell = 0.05*a0;

% use hexagonal lattice
[pos,ex,ey] = init_cellpos_hex(gridsize,gridsize);
dist = dist_mat(pos,gridsize,gridsize,ex,ey);
dist_vec = dist(1,:);

Son = 40;
K = 20;

dist_vec = dist_vec*a0;
r = dist_vec(dist_vec>0); % exclude self influence
fN = sum(Rcell*sum(exp(Rcell-r)./r)); % calculate signaling strength
gN = sum(sum((Rcell^2)*exp(2*(Rcell-r))./(r.^2))); % calculate noise variance strength

zoff = zeros(N+1,1);
zon = zeros(N+1,1);
factor = zeros(N+1,1);
for k = 0:N
    p = k/N;
    muON = Son + fN*(p*Son+(1-p));
    muOFF = 1 + fN*(p*Son+(1-p));
    sigmap = sqrt(p*(1-p)*gN)*(Son-1);
    zon(k+1) = (muON-K)/sigmap;
    zoff(k+1) = (K-muOFF)/sigmap;
    factor(k+1) = p*normcdf(zon(k+1)) + (1-p)*normcdf(zoff(k+1));
end

koff = N/(Son-1)*((K-1)/fN - 1)
kon = N/(Son-1)*((K-Son)/fN - 1)

figure(1)
plot(0:N,zoff,0:N,zon)

figure(2)
k = 0:N;
hold on
plot(0:N, (k'/N).^N.*factorial(N)./factorial(N-k')./factorial(k'))
plot(0:N, factor.^N.*factorial(N)./factorial(N-k')./factorial(k'))
hold off

figure(3)
hold on
plot(k,factor.^N)
plot(k,(k'/N).^N)
hold off

