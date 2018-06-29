%% Try to back-infer specific trajectory from final state of a trajectory
clear variables
close all

% load data
%{
fname = fullfile(pwd, 'figures\generate_pattern', 'test_trajectory.mat');
load(fname, 'cells_hist', 'a0', 'Rcell', 'Con', 'K', 't', 'gridsize');
N = gridsize^2;

[pos,ex,ey] = init_cellpos_hex_reverse(gridsize,gridsize);
dist = dist_mat(pos,gridsize,gridsize,ex,ey);
dist_vec = a0*dist(1,:);
r = dist_vec(dist_vec>0); % exclude self influence
fN = sum(sinh(Rcell)*sum(exp(Rcell-r)./r)); % calculate signaling strength
gN = sum(sum((sinh(Rcell)*exp(Rcell-r)./r).^2)); % calculate signaling strength
cell_type = zeros(N,1);

hin=figure();
for i=1:numel(cells_hist)
    k=waitforbuttonpress;
    update_cell_figure(hin, pos, 1, cells_hist{i}, cell_type, i-1);
end
%}
%% Run trajectory 
gridsize = 11;
N = gridsize^2;
a0 = 1.5;
Rcell = 0.2*a0;
Con = 14;
K = 16;
%}

%% Check trajectory agrees with predicted
cells_in = cells_hist{1};
hin=figure();
update_cell_figure(hin, pos, 1, cells_in, cell_type, i-1);
[cells, changed] = update_cells(cells_in, dist, Con, K, a0, Rcell);
while changed
    [cells, changed] = update_cells(cells, dist, Con, K, a0, Rcell);
    k=waitforbuttonpress;
    update_cell_figure(hin, pos, 1, cells, cell_type, i-1);
end
%% get Non at each time step
Non = zeros(t+1,1);
for i=1:t+1
    Non(i) = sum(cells_hist{i});
end
delNon = abs(Non(2:end)-Non(1:end-1));

figure();
plot(1:t, delNon);
%%
cells_pat = cells_hist{end};

% iteratively calculate preimage for a specific found cell
cells_past = {cells_pat};
nflips = 10000;
maxflips = 10;

finished = 0;
time = 0;
while ~finished
    cells_in = cells_past{end};
    found = 0;
    for flip=delNon(end-time)
        [cells_flipped, cells_accepted] = ...
            spin_flip_multiple__input_pattern(cells_in, dist, Con, K, a0, Rcell, flip, nflips);
        if sum(cells_accepted)>0
            idx = find(cells_accepted>0, 1); % store first found cell
            cells_past{end+1} = cells_flipped(idx);
            found = 1;
            time = time+1;
            break
        end
        fprintf('flips: %d \n', flip);
    end
    % if no pre-image found, stop
    if ~found
        disp('No pre-image found, finished.');
        finished = 1;
    end
end
