% Time evolution of a system without noise and without visualization

close all
clear all
warning off

data_path = '~/Dropbox/Matlab codes/data_onecelltype_entropy/no_self_comm';

% Parameters of the system
gridsize = 15;
N = gridsize^2;
a0 = 0.5;
Rcell = 0.2*a0;
p = 0.5;
iniON = round(p*N);
to_organize = 0;
n_run = 1;
% parameters
Son = 11;
K = 14;

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
    aux = round(to_organize*iniON);
    cells = zeros(N,1);
    cells(randperm(N,iniON-aux)) = 1;
    if aux > 0
        cells(find(cells==0,aux)) = 1;
    end
    load('tmp_data.mat','cells')
    % initialize figure
    %hin = figure(1);
    t = 0;
    I = [];
    Non = [];
    mom = [];
    cells_hist = {};
    %update_cell_figure(hin, pos, a0, cells, cell_type, t);
    cells_hist{end+1} = cells;
    Non(end+1) = sum(cells);
    I(end+1) = moranI(cells, a0*dist);
    [cells_out, changed, mom(end+1)] = update_cells_noselfcomm(cells, dist, Son, K, a0, Rcell);
    while changed && t < 200
        t = t+1;
        %k = waitforbuttonpress;
        %update_cell_figure(hin, pos, a0, cells_out, cell_type, t);
        cells_hist{end+1} = cells_out;
        I(end+1) = moranI(cells_out, a0*dist);
        Non(end+1) = sum(cells_out);
        cells = cells_out;
        [cells_out, changed, mom(end+1)] = update_cells_noselfcomm(cells, dist, Son, K, a0, Rcell);
    end
    
    fname_str = sprintf('N%d_n%d_neq_%d_a0%d_K_%d_Son_%d_t_%d', ...
        N, iniON, Non(end), 10*a0, K, Son, t);
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
