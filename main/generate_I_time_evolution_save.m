% Time evolution of a system without noise and without visualization
close all
clear all
warning off

% Parameters of the system
gridsize = 11;
N = gridsize^2;
a0 = 1.5;
Rcell = 0.2*a0;
K = 16;
Con = 50;

p0 = 0.5;
iniON = round(p0*N);
I0 = 0.2;
dI = 0.01;

n_run = 10;

% use hexagonal lattice
[dist, pos] = init_dist_hex(gridsize, gridsize);
dist_vec = a0*dist(1,:);
r = dist_vec(dist_vec>0); % exclude self influence
fN = sum(sinh(Rcell)*sum(exp(Rcell-r)./r)); % calculate signaling strength

% generate cell_type (0 case type 1, 1 case type 2)
cell_type = zeros(N,1); % all the same here

for tests = 1:n_run
    disp(tests)
    % initialize ON cells
    cells = zeros(N,1);
    cells(randperm(N,iniON)) = 1;
    
    test = false;
    counter = 0;
    while ~test && counter < 1000
        [cells_out, I_out, ~, test] = generate_I_new(cells, I0, I0+dI, dist, a0);
        counter = counter+1;
    end
    if test == false % if still not a good lattice
        disp('rejected');
        continue;
    end
    cells = cells_out;
    
    % initialize figure
    hin = figure(1);
    t = 0;
    I = [];
    Non = [];
    H = [];
    cells_hist = {};
    %corr = {}; % saves correlation function at different times
    
    update_cell_figure_withI(hin, pos, dist, a0, cells_out, cell_type, t);
    k = waitforbuttonpress;
    % store initial values
    Non(end+1) = sum(cells);
    I(end+1) = moranI(cells, a0*dist);
    cells_hist{end+1} = cells;
    %[corr{end+1}, ~] = correlation_func(dist, cells);
    [cells_out, changed, H(end+1)] = update_cells(cells, dist, Con, K, a0, Rcell);
    
    while changed
        t = t+1;
        %k = waitforbuttonpress;
        update_cell_figure_withI(hin, pos, dist, a0, cells_out, cell_type, t);
        cells_hist{end+1} = cells_out;
        %[corr{end+1}, ~] = correlation_func(dist, cells_out);
        I(end+1) = moranI(cells_out, a0*dist);
        Non(end+1) = sum(cells_out);
        cells = cells_out;
        [cells_out, changed, H(end+1)] = update_cells(cells, dist, Con, K, a0, Rcell);
    end
    
    fname_str = sprintf('N%d_n%d_neq_%d_a0_%d_K_%d_Con_%d_t_%d', ...
        N, iniON, Non(end), 10*a0, K, Con, t);
    i = 1;
    %fname = fullfile(pwd, 'data','dynamics', 'no_noise', ...
    %    strcat(fname_str,'-v',int2str(i),'.mat'));
    fname = fullfile(pwd, 'data', 'dynamics', 'generate_I',...
        strcat(fname_str,'-v',int2str(i),'.mat'));
    while exist(fname, 'file') == 2
        i=i+1;
        %fname = fullfile(pwd, 'data', 'dynamics', 'no_noise', ...
        %    strcat(fname_str,'-v',int2str(i),'.mat'));
        fname = fullfile(pwd, 'data', 'dynamics', 'generate_I',...
            strcat(fname_str,'-v',int2str(i),'.mat'));
    end
    save(fname)
end
