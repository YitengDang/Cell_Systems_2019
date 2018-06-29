% Generate trajectories of similar initial conditions by switching single
% cells 

% Time evolution of a system without noise and without visualization
close all
clear all
warning off

% Parameters of the system
gridsize = 11;
N = gridsize^2;
a0 = 1.5;
Rcell = 0.2*a0;
p = 0.4;
iniON = round(p*N);
K = 6;
Con = 15;

% use hexagonal lattice
[dist, pos] = init_dist_hex(gridsize, gridsize);
dist_vec = a0*dist(1,:);
r = dist_vec(dist_vec>0); % exclude self influence
fN = sum(sinh(Rcell)*sum(exp(Rcell-r)./r)); % calculate signaling strength

% generate cell_type (0 case type 1, 1 case type 2)
cell_type = zeros(N,1); % all the same here

% initialize ON cells
cells_ini = zeros(N,1);
cells_ini(randperm(N,iniON)) = 1;

for tests = 0:N
    disp(tests)
    % flip all configurations except first one
    cells = cells_ini;
    if tests~=0
       cells(tests) = ~cells(tests);
    end
    % Variables
    %hin = figure(1);
    t = 0;
    I = [];
    Non = [];
    mom = [];
    cells_hist = {};
    
    % store initial values
    %update_cell_figure(hin, pos, a0, cells, cell_type, t);
    %k = waitforbuttonpress;
    Non(end+1) = sum(cells);
    I(end+1) = moranI(cells, a0*dist);
    cells_hist{end+1} = cells;
    [cells_out, changed, mom(end+1)] = update_cells(cells, dist, Con, K, a0, Rcell);
    while changed
        t = t+1;
        %k = waitforbuttonpress;
        %update_cell_figure(hin, pos, a0, cells_out, cell_type, t);
        cells_hist{end+1} = cells_out;
        I(end+1) = moranI(cells_out, a0*dist);
        Non(end+1) = sum(cells_out);
        cells = cells_out;
        [cells_out, changed, mom(end+1)] = update_cells(cells, dist, Con, K, a0, Rcell);
    end
    
    % save trajectory
    qsave = 1;
    if qsave
        fname_str = sprintf('spin_flip_N%d_n%d_flip%d_neq_%d_a0%d_K_%d_Son_%d_t_%d', ...
            N, iniON, tests, Non(end), 10*a0, K, Con, t);
        i = 1;
        %fname = fullfile(pwd, 'data','dynamics', 'no_noise', ...
        %    strcat(fname_str,'-v',int2str(i),'.mat'));
        fname = fullfile(pwd, 'data', 'dynamics', 'single_spin_flip',...
            strcat(fname_str,'-v',int2str(i),'.mat'));
        while exist(fname, 'file') == 2
            i=i+1;
            %fname = fullfile(pwd, 'data', 'dynamics', 'no_noise', ...
            %    strcat(fname_str,'-v',int2str(i),'.mat'));
            fname = fullfile(pwd, 'data', 'dynamics', 'single_spin_flip',...
                strcat(fname_str,'-v',int2str(i),'.mat'));
        end
        save(fname)
    end
end
