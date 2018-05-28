% Time evolution of a system without noise and without visualization
% loop over different values of n (iniON)
close all
clear all
warning off

% Parameters of the system
gridsize = 15;
N = gridsize^2;
a0 = 0.5;
Rcell = 0.2*a0;
nrange = 0:N;
Klist = [18]; %[20 6 12]; %[3 6 10 15 17 20]; % (a0=0.5) [10 10 14 19 16 18];
Conlist = [11]; %[15 15 15]; %[24 21 21 20 14 14]; % (a0=0.5) [5 21 16 14 8 6];
I0_target = 0; % target initial I

% use hexagonal lattice
[dist, pos] = init_dist_hex(gridsize, gridsize);
dist_vec = a0*dist(1,:);
r = dist_vec(dist_vec>0); % exclude self influence
fN = sum(sinh(Rcell)*sum(exp(Rcell-r)./r)); % calculate signaling strength

% generate cell_type (0 case type 1, 1 case type 2)
cell_type = zeros(N,1); % all the same here
to_organize = 0;
%% Plot phenotype map
%{
Con_vec = 1:0.1:30;
K_vec = 1:0.1:20;
[out, map] = get_phenotype_map(a0, dist, Rcell, Con_vec, K_vec);

figure();
him=imagesc(K_vec, Con_vec, out);
set(him, 'AlphaData', out > 0); % invisible if not from any region
colormap(map);
set(gca, 'YDir', 'Normal');
xlabel('K');
ylabel('C_{ON}');
%}
%%
for i=1:numel(Klist)
    K = Klist(i);
    Con = Conlist(i);
for iniON = nrange
    disp(iniON)

    % initialize ON cells
    cells = zeros(N,1);
    cells(randperm(N,iniON)) = 1;
    dI = 0.01; % max. target offset
    [cells, I_ini, ~, test] = generate_I_new(cells, I0_target, I0_target+dI, dist, a0);
    
    % initialize figure
    %hin = figure(1);
    t = 0;
    I = [];
    Non = [];
    mom = [];
    cells_hist = {};
    %corr = {}; % saves correlation function at different times
    
    % store initial values
    Non(end+1) = sum(cells);
    I(end+1) = moranI(cells, a0*dist);
    cells_hist{end+1} = cells;
    %[corr{end+1}, ~] = correlation_func(dist, cells);
    [cells_out, changed, mom(end+1)] = update_cells(cells, dist, Con, K, a0, Rcell);
    while changed
        t = t+1;
        %k = waitforbuttonpress;
        %update_cell_figure(hin, pos, a0, cells_out, cell_type, t);
        cells_hist{end+1} = cells_out;
        %[corr{end+1}, ~] = correlation_func(dist, cells_out);
        I(end+1) = moranI(cells_out, a0*dist);
        Non(end+1) = sum(cells_out);
        cells = cells_out;
        [cells_out, changed, mom(end+1)] = update_cells(cells, dist, Con, K, a0, Rcell);
    end
    
    fname_str = sprintf('N%d_n%d_neq_%d_a0_%d_K_%d_Son_%d_t_%d', ...
        N, iniON, Non(end), 10*a0, K, Con, t);
    i = 1;
    fname = fullfile(pwd, 'data', 'trajectories', strcat(fname_str,'-v',int2str(i),'.mat'));
    while exist(fname, 'file') == 2
        i=i+1;
        fname = fullfile(pwd, 'data', 'trajectories', strcat(fname_str,'-v',int2str(i),'.mat'));
    end
    save(fname)
end
end