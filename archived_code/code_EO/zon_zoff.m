close all
clear all

% Plot the number of cells for which zon and zoff are within a range. zon
% and zoff are the variables inside the error function for calculation of
% the probabilities P(on->on) and P(off->off)

% Parameters of the system
gridsize = 11;
N = gridsize^2;
a0 = 0.5;
Rcell = 0.2*a0;

% use hexagonal lattice
[pos,ex,ey] = init_cellpos_hex(gridsize,gridsize);
dist = dist_mat(pos,gridsize,gridsize,ex,ey);

K_vec = 1:0.1:20;
Son_vec = 1:0.1:60;

[K, Son] = meshgrid(K_vec, Son_vec);

dist_vec = dist(1,:);
r = dist_vec(dist_vec>0); % exclude self influence
fN = sum(Rcell*sum(exp(Rcell-r)./r)); % calculate signaling strength
gN = sum(Rcell^2*sum(exp(2*(Rcell-r))./r.^2)); % calculate signaling strength

% value of z for which phi becomes relevant
l = 2;

av_n = @(x,S) fN*(x.*S + (1-x));
sigma = @(x,S) sqrt(x.*(1-x)*gN).*(S-1);

% zoff
a = (Son-1).^2.*(l^2*gN + fN^2);
b = (Son-1).*(2*fN^2-gN*l^2.*(Son-1)-2*(K-1)*fN);
c = (fN - (K-1)).^2;
delta = b.^2-4.*a.*c;
out = zeros(size(delta));
idx = delta>=0;
out(idx) = delta(idx);

plot_heat_map(K_vec, Son_vec, out, 'zoff delta', ...
    'K_T', 'S_{ON}', '\Delta_{OFF}')

p1 = zeros(size(delta));
p1(idx) = 0.5*(-b(idx)+sqrt(delta(idx)))./a(idx);
zoff = zeros(size(delta));
zoff(idx) = (K(idx) - 1 - av_n(p1(idx),Son(idx)))./sigma(p1(idx),Son(idx));

plot_heat_map(K_vec, Son_vec, zoff, 'zoff when p plus', ...
    'K_T', 'S_{ON}', 'z_{OFF}^{+}')

plot_heat_map(K_vec, Son_vec, p1, 'p plus zoff', ...
    'K_T', 'S_{ON}', '(k/N)_{+}')

p2 = zeros(size(delta));
p2(idx) = 0.5*(-b(idx)-sqrt(delta(idx)))./a(idx);
zoff = zeros(size(delta));
zoff(idx) = (K(idx) - 1 - av_n(p2(idx),Son(idx)))./sigma(p2(idx),Son(idx));

plot_heat_map(K_vec, Son_vec, zoff, 'zoff when p minus', ...
    'K_T', 'S_{ON}', 'z_{OFF}^{-}')

plot_heat_map(K_vec, Son_vec, p2, 'p minus zoff', ...
    'K_T', 'S_{ON}', '(k/N)_{-}')


% zon
a = (Son-1).^2.*(l^2*gN + fN^2);
b = (Son-1).*(2*fN^2-gN*l^2.*(Son-1)+2*(Son-K)*fN);
c = (fN + (Son-K)).^2;
delta = b.^2-4.*a.*c;
out = zeros(size(delta));
out(delta>=0) = delta(delta>=0);
figure('Name', 'zon delta')
colormap('hot');
imagesc(K_vec, Son_vec, out);
colorbar;
xlabel('K_T')
ylabel('S_{ON}')
