%% Try to back-infer specific trajectory from final state of a trajectory
clear variables
close all

%% Run trajectory 
gridsize = 15;
N = gridsize^2;
a0 = 1.5;
Rcell = 0.2*a0;
Con = 25;
K = 10;

[pos,ex,ey] = init_cellpos_hex_reverse(gridsize,gridsize);
dist = dist_mat(pos,gridsize,gridsize,ex,ey);
dist_vec = a0*dist(1,:);
r = dist_vec(dist_vec>0); % exclude self influence
fN = sum(sinh(Rcell)*sum(exp(Rcell-r)./r)); % calculate signaling strength
gN = sum(sum((sinh(Rcell)*exp(Rcell-r)./r).^2)); % calculate signaling strength
cell_type = zeros(N,1);

p0 = 0.5;
iniON = round(p0*N);
cells = zeros(N, 1);
cells(randperm(N, iniON)) = 1;

cells_hist = {};
hin=figure();
t=0;
update_cell_figure(hin, pos, 1, cells, cell_type, t);
[cells, changed] = update_cells(cells, dist, Con, K, a0, Rcell);
cells_hist{end+1} = cells;
while changed
    t = t+1;
    [cells, changed] = update_cells(cells, dist, Con, K, a0, Rcell);
    cells_hist{end+1} = cells;
    %k=waitforbuttonpress;
    update_cell_figure(hin, pos, 1, cells, cell_type, t);
end

%% get Non at each time step
Non = zeros(t+1,1);
for i=1:t+1
    Non(i) = sum(cells_hist{i});
end
delNon = abs(Non(2:end)-Non(1:end-1));

figure();
plot(1:t, delNon);

%% iteratively calculate preimage for a specific found cell
cells_pat = cells_hist{end}; % pattern to generate
cells_past = {cells_pat};
nflips = 100000;
maxflips = 10;

finished = 0;
rtime = 0; % reverse time, time = t - rtime
cells_rev = cells_hist{end};
while ~finished
    found = 0;
    for flip=delNon(end-rtime)
        disp('checkpoint 1');
        [cells_flipped, cells_accepted] = ...
            spin_flip_multiple__input_pattern(cells_rev, dist, Con, K, a0, Rcell, flip, nflips);
        if sum(cells_accepted)>0
            disp('checkpoint 2');
            idx = find(cells_accepted>0, 1); % store first found cell
            cells_temp = cells_flipped(idx);
            cells_past{end+1} = cells_temp{1};
            cells_rev = cells_temp{1};
            found = 1;
            rtime = rtime+1;
            break
        else 
            disp('No pre-image found, finished.');
            finished = 1;
        end
        fprintf('flips: %d \n', flip);
    end
    % if no pre-image found, stop
    if rtime == t
        disp('Time = 0 reached, finished.');
        finished = 1;
    end
end
%% show found trajectory
hin=figure();
for rtime=1:numel(cells_past)
    update_cell_figure(hin, pos, 1, cells_past{end-rtime+1}, cell_type, rtime-1);
    pause(1);
end