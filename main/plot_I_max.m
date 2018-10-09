clear variables
close all
%%
gz = 15;
N = gz^2;
a0 = 1.5;
Rcell = 0.2*a0;

% Calculate fN
[dist, pos] = init_dist_hex(gz, gz);
dist_vec = a0*dist(1,:);
r = dist_vec(dist_vec>0); % exclude self influence
fN = sum(sinh(Rcell)*sum(exp(Rcell-r)./r)); % calculate signaling strength

pv = (0:N)/N;
figure;
hold on
fa0 = sinh(Rcell)*exp(Rcell-a0)/a0;
Imax = @(p) (6 - 4./sqrt(p*N) - 6./(p*N) - 6*p)*fa0./(1-p)/fN;
Imax2 = @(p) (6*p - 4./sqrt((1-p)*N) - 6./((1-p)*N))*fa0./p/fN;
plot(pv, Imax(pv), 'b')
plot(pv, max(Imax(pv), Imax2(pv)), '--b', 'Linewidth', 1.5)
ylim([-0.05 1])
xlabel('p');
ylabel('I');