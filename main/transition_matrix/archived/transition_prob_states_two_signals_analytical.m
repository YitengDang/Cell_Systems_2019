% Calculates the state transition probabilities for two signals through
% analytical theory
% pij: for transitions between (i,j) states
clear variables
close all
warning off
set(0, 'defaulttextinterpreter', 'latex');

% Parameters of the system
gz = 8;
N = gz^2;
a0 = 1.5;
rcell = 0.2;
Rcell = rcell*a0;

M_int = [1 -1; 1 0];
Con = [35 15];
Coff = [1 1];
K = [8 25; 35 0];
lambda = [1 1.2];
hill = Inf;
noise = 0;

% (1) original hexagonal lattice
%[pos,ex,ey] = init_cellpos_hex(gridsize,gridsize);
%dist = dist_mat(pos,gridsize,gridsize,ex,ey);
% (2) new lattice
mcsteps = 0;
nodisplay = 1;
[pos, dist] = initial_cells_random_markov_periodic(gz, mcsteps, rcell, nodisplay);

% Calculate interaction strength
dist_vec = a0*dist(1,:);
r = dist_vec(dist_vec>0); % exclude self influence
fN = zeros(1, 2);
gN = zeros(1, 2);
for i=1:2
    fN(i) = sum(sinh(Rcell) * sum(exp((Rcell-r)./lambda(i)).*(lambda(i)./r)) ); % calculate signaling strength
    gN(i) = sum(sinh(Rcell)^2 * sum(exp(2*(Rcell-r)./lambda(i)).*(lambda(i)./r).^2 ) ); % calculate noise variance strength
end

n_trials = 10; % simulations to consider

% Initial state
iniON = zeros(2);
iniON(1,1) = 15;
iniON(1,2) = 12;
iniON(2,1) = 18;
iniON(2,2) = N - iniON(1,1) - iniON(1,2) - iniON(2,1); %16;

if sum(sum(iniON))~=N
    warning('Wrong probabilities');
end
%% Calculate transition matrix
W = zeros(4, 4); %transition matrix
for trial=1:n_trials
    % generate initial lattice

    % subdivide indices
    idx = cell(3,1);
    idx_cells_not00 = randperm(N, N-iniON(1,1)); % cell indices of all cells not (0,0)
    part_idx = [1 iniON(1,2) iniON(1,2)+iniON(2,1) iniON(1,2)+iniON(2,1)+iniON(2,2)]; % partition indices
    idx{1} = idx_cells_not00(part_idx(1):part_idx(2));
    idx{2} = idx_cells_not00(part_idx(2)+1:part_idx(3));
    idx{3} = idx_cells_not00(part_idx(3)+1:part_idx(4));

    % Check that sets have the right number of elems
    check = (numel(idx{1})==iniON(1,2));
    check = check*(numel(idx{2})==iniON(2,1));
    check = check*(numel(idx{3})==iniON(2,2));
    % Check that sets have disjoint elements
    check = check*(numel(union(union(idx{1}, idx{2}), idx{3}))==N-iniON(1,1)); 
    fprintf('Check passed? %d \n', check);

    % turn genes ON
    cellsIn = zeros(N, 2);
    cellsIn(idx{1}, 2) = 1;
    cellsIn(idx{2}, 1) = 1;
    cellsIn(idx{3}, 1) = 1;
    cellsIn(idx{3}, 2) = 1;

    % Update cells
    t = 0;
    [cellsOut, changed] = update_cells_two_signals_multiply_finite_Hill(cellsIn, dist, M_int, a0,...
            Rcell, Con, Coff, K, lambda, hill, noise);
        
    cellsIn_idx = cellsIn*[2; 1];
    cellsOut_idx = cellsOut*[2; 1];
    for i=1:N
        W(cellsOut_idx(i)+1, cellsIn_idx(i)+1) = W(cellsOut_idx(i)+1, cellsIn_idx(i)+1) + 1;
    end
end
W = W/(n_trials*N);

%% Plot matrix
h = figure(1);
imagesc(1:4, 1:4, W)
c = colorbar;
set(gca, 'Ydir', 'normal', 'FontSize', 20)
xlabel('$$p_1(t+1)$$', 'FontSize', 24)
ylabel('$$p_2(t+1)$$', 'FontSize', 24)
c.Label.String = 'Probability';
%}
