% Time evolution of a system with noise and without visualization

close all
clear all
warning off

% Parameters of the system
gridsize = 15;
N = gridsize^2;
a0 = 1.5;
Rcell = 0.2*a0;
p = 0.4;
iniON = round(p*N);
to_organize = 0;
n_run = 100;
noise = 5;
% parameters
Son = 21.7;
K = 15.8;
%tmax = 100;

data_path = 'D:\eduardopavinat\Dropbox\Matlab codes\data_onecelltype_entropy';

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
    [cells_out, changed, mom(end+1)] = update_cells_noise(cells, dist, Son, K, a0, Rcell, noise);
    % run until stabilize
    cont = true;
    while cont
        t = t+1;
        %k = waitforbuttonpress;
        %update_cell_figure(hin, pos, a0, cells_out, cell_type, t);
        cells_hist{end+1} = cells_out;
        I(end+1) = moranI(cells_out, a0*dist);
        Non(end+1) = sum(cells_out);
        cells = cells_out;
        [cells_out, changed, mom(end+1)] = update_cells_noise(cells, dist, Son, K, a0, Rcell, noise);
        if t >= 10
            if std(Non(end-9:end)) <= 1 && std(I(end-9:end)) <=0.01
                cont = false;
            end
        end
    end
    
    fname_str = sprintf('N%d_n%d_neq_%d_a0%d_K_%d_Son_%d_noise_%d_t_%d', ...
        N, iniON, Non(end), round(10*a0), round(K), round(Son), round(10*noise), t);
    i = 1;
    % check if filename exists, if not change the end of the filename
    fname = fullfile(data_path, 'dynamics_noise', ...
        strcat(fname_str,'-v',int2str(i),'.mat'));
    while exist(fname, 'file') == 2
        i=i+1;
        fname = fullfile(data_path,'dynamics_noise', ...
            strcat(fname_str,'-v',int2str(i),'.mat'));
    end
    
    save(fname)
end
