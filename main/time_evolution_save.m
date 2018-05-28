% Time evolution of a system without noise and without visualization
close all
clear all
warning off
set(0, 'defaulttextinterpreter', 'latex');

% Parameters of the system
gridsize = 15;
N = gridsize^2;
a0 = 0.5;
Rcell = 0.2*a0;

% circuit parameters
Con = 11;
K = 18;

% initial conditions
p0_all = (0:5:N)/N;
%p0 = 0.5;
%iniON = round(p0*N);
I0 = 0;
to_organize = 0;
%n_run = 50;

% use hexagonal lattice
[dist, pos] = init_dist_hex(gridsize, gridsize);
dist_vec = a0*dist(1,:);
r = dist_vec(dist_vec>0); % exclude self influence
fN = sum(sinh(Rcell)*sum(exp(Rcell-r)./r)); % calculate signaling strength

% generate cell_type (0 case type 1, 1 case type 2)
cell_type = zeros(N,1); % all the same here

for tests = 1:numel(p0_all) %n_run
    disp(tests);
    p0 = p0_all(tests);
    iniON = round(p0*N);

    % initialize ON cells
    aux = round(to_organize*iniON);
    cells = zeros(N,1);
    cells(randperm(N,iniON-aux)) = 1;
    if aux > 0
        cells(find(cells==0,aux)) = 1;
    end
    [cells] = generate_I_new(cells, I0, I0+0.01, dist, a0);
  
    % initialize figure
    %hin = figure(1);
    t = 0;
    I = [];
    Non = [];
    h = [];
    cells_hist = {};
    %corr = {}; % saves correlation function at different times
    
    % store initial values
    Non(end+1) = sum(cells);
    I(end+1) = moranI(cells, a0*dist);
    cells_hist{end+1} = cells;
    %[corr{end+1}, ~] = correlation_func(dist, cells);
    [cells_out, changed, h(end+1)] = update_cells(cells, dist, Con, K, a0, Rcell);
    while changed
        t = t+1;
        %k = waitforbuttonpress;
        %update_cell_figure(hin, pos, a0, cells_out, cell_type, t);
        cells_hist{end+1} = cells_out;
        %[corr{end+1}, ~] = correlation_func(dist, cells_out);
        I(end+1) = moranI(cells_out, a0*dist);
        Non(end+1) = sum(cells_out);
        cells = cells_out;
        [cells_out, changed, h(end+1)] = update_cells(cells, dist, Con, K, a0, Rcell);
    end
    
    fname_str = strrep(sprintf('N%d_n%d_I%.1f_neq_%d_a0_%.2f_K_%d_Con_%d_t_%d', ...
        N, iniON, I0, Non(end), a0, K, Con, t), '.', 'p');
    i = 1;
    %fname = fullfile(pwd, 'data','dynamics', 'no_noise', ...
    %    strcat(fname_str,'-v',int2str(i),'.mat'));
    fname = fullfile(pwd, 'data', 'dynamics', 'no_noise',...
         strcat(fname_str,'-v',int2str(i),'.mat'));
    while exist(fname, 'file') == 2
        i=i+1;
        %fname = fullfile(pwd, 'data', 'dynamics', 'no_noise', ...
        %    strcat(fname_str,'-v',int2str(i),'.mat'));
        fname = fullfile(pwd, 'data', 'dynamics', 'no_noise',...
            strcat(fname_str,'-v',int2str(i),'.mat'));
    end
    
    save(fname)
end
