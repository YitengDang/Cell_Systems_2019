% Generate trajectories of similar initial conditions by switching single
% cells 

% Time evolution of a system without noise and without visualization
close all
clear all
warning off

% Parameters
gridsize = 11;
N = gridsize^2;
a0 = 5.4;
Rcell = 0.2*a0;
p = 0.1;
iniON = round(p*N);
K = 8;
Con = 16;
initialID = 'binary';
hill = 2;
prec = 8;

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

for flip = 0:N
    disp(flip)
    % flip all configurations except first one
    cells = cells_ini;
    if flip~=0
       cells(flip) = ~cells(flip);
    end
    % Variables
    %hin = figure(1);
    t = 0;
    I2 = [];
    Theta = [];
    Non = [];
    H = [];
    cells_hist = {};
    
    % store initial values
    %update_cell_figure(hin, pos, a0, cells, cell_type, t);
    %k = waitforbuttonpress;
    Non(end+1) = sum(cells);
    [I2(end+1), Theta(end+1)] = moranI_new(cells, a0*dist);
    cells_hist{end+1} = cells;
    [cells_out, changed, H(end+1)] = update_cells_continuum(cells, dist, Con, K, a0, Rcell, hill, prec);
    while changed
        t = t+1;
        %k = waitforbuttonpress;
        %update_cell_figure(hin, pos, a0, cells_out, cell_type, t);
        cells_hist{end+1} = cells_out;
        [I2(end+1), Theta(end+1)] = moranI_new(cells_out, a0*dist);
        Non(end+1) = sum(cells_out);
        cells = cells_out;
        [cells_out, changed, H(end+1)] = update_cells_continuum(cells, dist, Con, K, a0, Rcell, hill, prec);
    end
    
    % save trajectory
    qsave = 1;
    if qsave
        fname_str = strrep(sprintf('perturb_ini_N%d_pini%.3f_flip%d_a0_%.1f_K_%d_Con_%d_hill_%.2f', ...
            N, p, flip, a0, K, Con, hill), '.', 'p');
        i = 1;
        %fname = fullfile(pwd, 'data','dynamics', 'no_noise', ...
        %    strcat(fname_str,'-v',int2str(i),'.mat'));
        fname = fullfile(pwd, 'data', 'perturb_ini',...
            strcat(fname_str,'-v',int2str(i),'.mat'));
        while exist(fname, 'file') == 2
            i=i+1;
            %fname = fullfile(pwd, 'data', 'dynamics', 'no_noise', ...
            %    strcat(fname_str,'-v',int2str(i),'.mat'));
            fname = fullfile(pwd, 'data', 'perturb_ini',...
                strcat(fname_str,'-v',int2str(i),'.mat'));
        end
        save(fname)
    end
end
