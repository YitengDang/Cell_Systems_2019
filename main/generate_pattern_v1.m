%% Design a desired pattern and try to find parameters for which it is stable
% v1: old
clear variables
close all
%% Input pattern
input = ...
    [
        0	0	0	0	0	0	0	0	0	0;
        0	1	0	1	0	1	0	0	0	1;
        0	1	0	1	0	0	1	0	1	0;
        0	1	0	1	0	0	0	1	0	0;
        0	1	1	1	0	0	0	1	0	0;
        0	1	0	1	0	0	0	1	0	0;
        0	1	0	1	0	0	0	1	0	0;
        0	1	0	1	0	0	0	1	0	0;
        0	1	0	1	0	0	0	1	0	0;
        0	0	0	0	0	0	0	0	0	0;
    ];

%% Variables
a0 = 0.5; % try different values
Rcell = 0.2*a0;

cells = input(:);
N = numel(cells);
if ~mod(N,1)
    gridsize = sqrt(N);
else 
    disp('Error');
end
cell_type = zeros(N, 1);

[pos,ex,ey] = init_cellpos_hex_reverse(gridsize,gridsize);
dist = dist_mat(pos,gridsize,gridsize,ex,ey);
dist_vec = a0*dist(1,:);
r = dist_vec(dist_vec>0); % exclude self influence
fN = sum(sinh(Rcell)*sum(exp(Rcell-r)./r)); % calculate signaling strength
gN = sum(sum((sinh(Rcell)*exp(Rcell-r)./r).^2)); % calculate signaling strength

% Plot pattern
hin=figure();
update_cell_figure(hin, pos, 1, cells, cell_type, 0)

%% Look for parameters for which pattern is stable
p = sum(cells)/N;
[I, Theta] = moranI(cells, a0*dist);

% Plot Peq map for fixed I, Theta against K, Con
Con_all = 1:40;
K_all = 1:40;
[Con_mesh, K_mesh] = meshgrid(Con_all, K_all);

muon = fN*(Con_mesh.*p + 1 - p + (Con_mesh-1).*(1-p).*I);
muoff = fN*(Con_mesh.*p + 1 - p - (Con_mesh-1).*p.*I);
kappa = sqrt((Con_mesh-1).^2*(gN)*p.*(1-p));
sigmaon = kappa;
sigmaoff = kappa;

zon = (K_mesh - Con_mesh - muon)./sigmaon;
zoff = (K_mesh - 1 - muoff)./sigmaoff;

Poffoff = normcdf(zoff);
Ponon = 1-normcdf(zon);

pe = (Ponon.^p .* Poffoff.^(1-p)).^N;

figure();
imagesc(Con_all, K_all, pe);
set(gca, 'YDir', 'normal', 'FontSize', 24);
c=colorbar;

%% Find K, Con where P_eq is maximal and test whether pattern is stable
[ix, iy] = find(pe == max(max(pe)));
fprintf('Max P_eq = %.3f \n', pe(ix, iy));
K = K_all(ix);
Con = Con_all(iy);

[cells_pat, changed, ~] = update_cells(cells, dist, Con, K, a0, Rcell);
if ~changed % stable configuration
    disp('stable pattern!');
end
%% algorithm flipping 1 cell at a time
test = false;
idx = 1;
while ~test && idx<N+1
    disp(idx);
    % flip random cell
    cells_flip = cells_pat;
    cells_flip(idx) = ~cells_flip(idx);

    hin=figure(1);
    update_cell_figure(hin, pos, 1, cells_flip, cell_type, 0)
    k=waitforbuttonpress;

    % check whether the system evolves back to config
    [cells_temp, changed, ~] = update_cells(cells_flip, dist, Con, K, a0, Rcell);
    if all(cells_temp==cells_pat)
        test = true;
    end
    idx = idx+1;
end
