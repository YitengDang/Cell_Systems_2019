%% Design a desired pattern and try to find parameters for which it is stable
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

cells_pat = input(:);
N = numel(cells_pat);
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
hin=figure(1);
update_cell_figure(hin, pos, 1, cells_pat, cell_type, 0)

%% Look for parameters for which pattern is stable
p = sum(cells_pat)/N;
[I, Theta] = moranI(cells_pat, a0*dist);

% Plot Peq map for fixed I, Theta against K, Con
Con_all = 1:30;
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
%% Find K, Con where P_eq is large and test whether pattern is stable
% get list of sorted P_eq
[pe_sorted, I] = sort(pe(:), 'descend');
idx = pe_sorted > 0;
[Kidx, Conidx] = ind2sub([numel(K_all) numel(Con_all)], I(idx));

% test for all K, Con that satisfy constraint
stability_list = zeros(numel(Kidx), 1);
for i=1:numel(Kidx)
    K = Kidx(i);
    Con = Conidx(i);
    [~, changed, ~] = update_cells(cells_pat, dist, Con, K, a0, Rcell);
    stability_list(i) = ~changed;
    if ~changed % stable configuration
        disp('stable pattern!');
    end
end
K_stable = Kidx(stability_list>0);
Con_stable = Conidx(stability_list>0);

h=figure(2);
hold on
imagesc(Con_all, K_all, pe);
% plot chosen points
plot(Conidx, Kidx, 'gx');
% plot stable points
%idx = stability_list>0;
%plot(Conidx(idx), Kidx(idx), 'ro');
plot(Con_stable, K_stable, 'ro');
% plot style
set(gca, 'YDir', 'normal', 'FontSize', 24);
c=colorbar;
xlabel('Con');
ylabel('K')
xlim([Con_all(1)-1 Con_all(end)+1]);
ylim([K_all(1)-1 K_all(end)+1]);
%% algorithm finding prior states by flipping spins
%{
% (1) for specific K, Con
K = K_stable(1);
Con = Con_stable(1);
nflips = 10000;
frac_accepted = zeros(N, 1);
for flip=1:10
    [cells_flipped, cells_accepted] = ...
        spin_flip_multiple__input_pattern(cells_pat, dist, Con, K, a0, Rcell, flip, nflips);
    frac_accepted(flip) = sum(cells_accepted)/numel(cells_accepted);
    fprintf('test %d \n', flip);
    fprintf('Fraction accepted: %.3f \n', frac_accepted(flip));
end

figure();
plot(1:N, frac_accepted);
xlabel('number of flips');
ylabel('fraction accepted');
%}
%% (2) for all stable K, Con
maxflips = 1;
frac_accepted = zeros(numel(K_stable), maxflips);
cells_preimage = cell(numel(K_stable), maxflips);
for i=1:numel(K_stable)
    K = K_stable(i);
    Con = Con_stable(i);
    fprintf('K = %d, Con = %d \n', K, Con);
    nflips = 10000;    
    for flip=1:maxflips
        [cells_flipped, cells_accepted] = ...
            spin_flip_multiple__input_pattern(cells_pat, dist, Con, K, a0, Rcell, flip, nflips);
        % store accepted cells
        idx = find(cells_accepted>0);
        cells_add = cells_flipped(idx);
        cells_preimage{i, flip} = cells_add;
        % update accepted fraction
        frac_accepted(i, flip) = sum(cells_accepted)/numel(cells_accepted);
        fprintf('test %d flips \n', flip);
        fprintf('Fraction accepted: %.3f \n', frac_accepted(flip));
    end
end

%% plot result
h3=figure(3);
imagesc(1:maxflips, 1:numel(K_stable), frac_accepted*nflips)
ylabel('Parameter value');
xlabel('Number of flips');
c=colorbar;
%% Get pre-images and associated parameters
% Get parameters for preimage
[param_out, flip_out] = find(frac_accepted>0);

% create cell array of pre-images with their respective parameters
cells_preimage2 = {};
params_preimage2 = [];
for i=1:numel(param_out)
    k = numel(cells_preimage{param_out(i), flip_out(i)} );
    disp(k);
    this_Con = Con_stable(param_out(i));
    this_K = K_stable(param_out(i));
    params_preimage2 = [params_preimage2; repmat([this_Con this_K], k, 1)];
    cells_preimage2 = [cells_preimage2; cells_preimage{param_out(i), flip_out(i)}];
end

%%
% Plot pre-image configurations
for idx=1:size(params_preimage2, 1)
disp(idx);
found = 0;
%idx = 2;
cells_preim_in = cells_preimage2{idx};
K = params_preimage2(idx, 2);
Con = params_preimage2(idx, 1);
cells_past = {cells_pat cells_preim_in};

% plot pre-image pattern and check that pattern indeed goes to correct result
hin=figure(4);
update_cell_figure(hin, pos, 1, cells_preim_in, cell_type, 0)
[cells, changed, ~] = update_cells(cells_preim_in, dist, Con, K, a0, Rcell);
k=waitforbuttonpress;
update_cell_figure(hin, pos, 1, cells, cell_type, 1)

% iteratively calculate preimage for a specific found cell
nflips = 10000;
maxflips = 10;
for flip=1:maxflips
    [cells_flipped, cells_accepted] = ...
        spin_flip_multiple__input_pattern(cells_preim_in, dist, Con, K, a0, Rcell, flip, nflips);
    if sum(cells_accepted)>0
        idx = find(cells_accepted>0, 1); % store first found cell
        cells_past{end+1} = cells_flipped(idx);
        found = 1;
        break
    end
    fprintf('flips: %d \n', flip);
end
if found
    break
end

end
