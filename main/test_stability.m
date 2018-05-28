% Calculate if a determined configuration is stable for different
% combinations of parameters

clear variables
close all

% Parameters
gridsize = 15;
a0 = 0.5;
Rcell = 0.2*a0;

[dist, pos] = init_dist_hex(gridsize, gridsize);

N = gridsize^2;
Con_vec = linspace(1,30,59);
K_vec = linspace(1,20,39);

M = sinh(Rcell)*exp(Rcell-a0*dist)./(a0*dist);
M(1:N+1:N^2) = 0;
aux = M*M';
gN = aux(1,1);
xiN = sum(aux(1,2:end));
fN = sum(M(1,:));

%%% --- Configuration --- %%%
cells = zeros(N, 1);

% stripe
% stripe_size = 3;
% cells(1:stripe_size*gridsize) = 1;

% circle
radius_size = 0.4;
centerx = gridsize/2*a0 + gridsize/4*a0;
centery = gridsize/4*a0*sqrt(3);

pos_mod = a0*pos;
pos_mod(:,1) = pos_mod(:,1) - centerx;
pos_mod(:,2) = pos_mod(:,2) - centery;

r = sqrt(sum(pos_mod.^2, 2));

cells(r<gridsize*a0*radius_size) = 1;

%%% --------------------- %%%

I = moranI(cells, a0*dist);
p = sum(cells)/N;
h1 = figure(1);
update_cell_figure_withI(h1, pos, dist, a0, cells, zeros(N,1), 0)

out = zeros(numel(K_vec), numel(Con_vec));
for i = 1:numel(K_vec)
    for j = 1:numel(Con_vec)
        [~, changed, mom] = ...
            update_cells(cells, dist, Con_vec(j), K_vec(i), a0, Rcell);
        if ~changed
            out(i, j) = 1;
        end
    end
end

[K, Con] = meshgrid(K_vec, Con_vec);
muon = fN*(Con*p + 1 - p + (Con-1)*(1-p)*I);
muoff = fN*(Con*p + 1 - p - (Con-1)*p*I);
sigma = sqrt((Con-1).^2*(gN)*p*(1-p));
zon = (K - Con - muon)./sigma;
zoff = (K - 1 - muoff)./sigma;
Poffoff = normcdf(zoff);
Ponon = 1-normcdf(zon);
pe = transpose((Ponon.^p.*Poffoff.^(1-p)).^N);

h2 = figure(2);
im_fig = imagesc(K_vec, Con_vec, out');
set(im_fig, 'AlphaData', out' ~= 0);
set(gca, 'YDir', 'normal', 'FontSize', 20)
hold on
tmp = find(pe > 1 - 1/N);
[I, J] = ind2sub(size(pe), tmp);
x = [];
y1 = [];
y2 = [];
for i = 1:numel(K_vec)
    tmp = J(I==i);
    if ~isempty(tmp)
        x(end+1) = K_vec(i);
        y1(end+1) = Con_vec(min(tmp));
        y2(end+1) = Con_vec(max(tmp));
    end
end
plot(x, y1, x, y2)
hold off


