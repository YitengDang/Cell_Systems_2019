% Two cell dynamics until equilibrium (no noise)
close all
clear all
warning off

for iii = 1:100
% Parameters of the system
gridsize = 20;
N = gridsize^2;
a0 = 0.5;
Rcell = 0.2*a0;
N1 = round(0.3*N);
N2 = N-N1;
Nv = [N1 N2];
p1 = 0.5;
p2 = 0.5;
n1 = round(p1*N1);
n2 = round(p2*N2);

% use hexagonal lattice
[dist, pos] = init_dist_hex(gridsize, gridsize);
dist_vec = dist(1,:);

K1 =16;
K2 = 7;
Son1 = 8;
Son2 = 5;

% Randomize the location of type 1 and 2 cells
idx1 = randperm(N,N1);
idx2 = setdiff(1:N, idx1);

% f_ij matrix M
M = zeros(N, N);
M(dist>0) = sinh(Rcell).*exp(Rcell-a0*dist(dist>0)) ...
    ./(a0*dist(dist>0));

f11 = sum(sum(M(idx1,idx1)))/N;
f12 = sum(sum(M(idx1,idx2)))/N;
f22 = sum(sum(M(idx2,idx2)))/N;

F = zeros(2,2);
F(1,1) = f11;
F(1,2) = f12;
F(2,1) = f12;
F(2,2) = f22;

% Set the genetic circuit parameters of the cell types
Kv = ones(N,1);
Sonv = ones(N,1);
Kv(idx1) = K1;
Kv(idx2) = K2;
Sonv(idx1) = Son1;
Sonv(idx2) = Son2;

% Calculate B's
S = [Son1 Son2];
K = [K1 K2];
B = 0.5*(S+1) - K + N./Nv .* (0.5*(S+1)*F);

% Set cells randomly ON
cells = zeros(N,1);
cells(idx1(randperm(N1,n1))) = 1;
cells(idx2(randperm(N2,n2))) = 1;
cells_hist = {};
cells_hist{end+1} = cells;
t = 0;
[cells, changed] = update_cells(cells, dist, Sonv, Kv, a0, Rcell);
while changed
    t = t+1;
    cells_hist{end+1} = cells;
    [cells, changed] = update_cells(cells, dist, Sonv, Kv, a0, Rcell);
end

fname_str = sprintf('N1-%d_n1-%d_N2-%d_n2-%d_a0-%d_K1-%d_K2-%d_Son1-%d_Son2-%d', ...
    N1, n1, N2, n2, round(10*a0), round(K1), round(K2), round(Son1), round(Son2));
i = 1;
fname = fullfile(data_path, 'dynamics_nonoise', ...
    strcat(fname_str,'-v',int2str(i),'.mat'));
while exist(fname, 'file') == 2
    i=i+1;
    fname = fullfile(data_path,'dynamics_nonoise', ...
        strcat(fname_str,'-v',int2str(i),'.mat'));
end

save(fname)
end