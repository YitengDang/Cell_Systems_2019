% Time evolution of a system without noise, with burst, without visualization

% Takes only trajectories that start in autonomy

% Results to be analyzed with 
% plot_saved_paths_with_hamiltonian_map_burst.m (Plot paths on H-map)
% energy_change_exact_burst.m (work-energy relationship)
close all
clear all
warning off

% Lattice parameters
gridsize = 11;
N = gridsize^2;
a0 = 0.5;
Rcell = 0.2*a0;

% initial conditions
p = 0.65;
iniON = round(p*N);

% circuit parameters
K = 16;
Son = 8;

% burst parameters
B = 200; %burst sizes
%bt = 1; %burst times, assume min(bt)>0

% simulation parameters
to_organize = 0;
n_run = 5;

% use hexagonal lattice
[dist, pos] = init_dist_hex(gridsize, gridsize);

dist_vec = a0*dist(1,:);
r = dist_vec(dist_vec>0); % exclude self influence
fN = sum(sinh(Rcell)*sum(exp(Rcell-r)./r)); % calculate signaling strength

% generate cell_type (0 case type 1, 1 case type 2)
cell_type = zeros(N,1); % all the same here
%%
for tests = 1:n_run
    disp(strcat('run ', int2str(tests)));
    % initialize ON cells
    aux = round(to_organize*iniON);
    cells = zeros(N,1);
    cells(randperm(N,iniON-aux)) = 1;
    if aux > 0
        cells(find(cells==0,aux)) = 1;
    end
    % initialize variables
    t = 0;
    I = [];
    Non = [];
    mom = [];
    cells_hist = {};
    % store initial configuration
    cells_hist{end+1} = cells;
    Non(end+1) = sum(cells);
    I(end+1) = moranI(cells, a0*dist);
    % initial step: check whether cells are in equilibrium
    [cells_out, changed, mom(end+1)] = update_cells(cells, dist, Son, K, a0, Rcell);
    if changed %take only configs that start in steady state
        disp('Initial state not stationary!');
        continue;
    else
        % if in equilibrium, update using burst
        [cells_out, changed, mom(end+1)] = update_cells_burst(cells, dist, Son, K, a0, Rcell, B);
    end
    while changed
        % currently not considering burst at t>0
        %{ 
        this_bt=find(bt==t, 1);
        this_B = 0;
        if ~isempty(this_bt)
            this_B = B(this_bt); 
        end
        disp(strcat('B=',num2str(this_B)));
        %} 
        t = t+1;
        cells_hist{end+1} = cells_out;
        I(end+1) = moranI(cells_out, a0*dist);
        Non(end+1) = sum(cells_out);
        cells = cells_out;
        [cells_out, changed, mom(end+1)] = update_cells(cells, dist, Son, K, a0, Rcell);
    end
    
    % save file
    fname_str = sprintf('N%d_n%d_neq%d_a0%d_K%d_Son%d_t%d_B%d', ...
        N, iniON, Non(end), 10*a0, K, Son, t, B);
    i = 1;
    fname = fullfile(pwd, 'data','dynamics', 'burst', ...
        strcat(fname_str,'-v',int2str(i),'.mat'));
    while exist(fname, 'file') == 2
        i=i+1;
        fname = fullfile(pwd, 'data', 'dynamics', 'burst', ...
            strcat(fname_str,'-v',int2str(i),'.mat'));
    end
    save(fname)
end
