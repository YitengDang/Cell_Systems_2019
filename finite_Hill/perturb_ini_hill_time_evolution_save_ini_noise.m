% Generate trajectories of similar initial conditions by switching single
% cells 

% Time evolution of a system without noise and without visualization
close all
clear variables
warning off

% Parameters
gridsize = 11;
N = gridsize^2;
a0 = 5.4;
Rcell = 0.2*a0;
p = 0.3;
iniON = round(p*N);
K = 8;
Con = 16;
hill = 2;
prec = 8;
%ini_noise = 0.1;
ini_noise_all = 10.^(-2:0.2:0);

initialID = 'normal_noise'; %'binary';
ntrials = 100;

save_folder_data = 'H:\My Documents\Multicellular automaton\figures\finite_Hill\perturb_ini\data_binary_constraint_v2\vs_perturbation_strength';

% use hexagonal lattice
[dist, pos] = init_dist_hex(gridsize, gridsize);
dist_vec = a0*dist(1,:);
r = dist_vec(dist_vec>0); % exclude self influence

fN = sum(sinh(Rcell)*sum(exp(Rcell-r)./r)); % calculate signaling strength
% Find uniform lattice fixed points
fp = zeros(3, 1);
x0 = [0.03 0.2 0.65]; %estimates based on previous graph
%hfunc = @update_function_uniform;
hfunc2 = @(x) ((1+fN)*((Con-1)*x + 1))^hill/(K^hill+ (((1+fN)*(Con-1)*x + 1))^hill) - x;
for idx=1:3
    fp(idx) = fzero(hfunc2, x0(idx));
end
    
% generate cell_type (0 case type 1, 1 case type 2)
cell_type = zeros(N,1); % all the same here

% loop over ini_noise
for i=1:numel(ini_noise_all)
    ini_noise = ini_noise_all(i);
    % initialize ON cells
    % -- binary cells -- 
    %{
    cells_ini = zeros(N,1);
    cells_ini(randperm(N,iniON)) = 1;
    %}
    % -- binary constraint v2 --
    cells_ini = init_random_cells_binary_constraint_v2(N, iniON, fp);

    for trial = 0:ntrials
        disp(trial)
        % flip all configurations except first one
        cells = cells_ini;
        if trial>0
            cells = generate_pattern_noise(cells, ini_noise); % add noise to original pattern
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
            fname_str = strrep(sprintf('perturb_ini_N%d_pini%.3f_a0_%.1f_K_%d_Con_%d_hill_%.2f_ini_noise_%.2f', ...
                N, p, a0, K, Con, hill, ini_noise), '.', 'p');
            v = 0;
            fname = fullfile(save_folder_data, strcat(fname_str,'-v',int2str(v),'.mat'));
            %fname = fullfile(pwd, 'data', 'perturb_ini',...
            %    strcat(fname_str,'-v',int2str(i),'.mat'));

            while exist(fname, 'file') == 2
                v=v+1;
                fname = fullfile(save_folder_data, strcat(fname_str,'-v',int2str(v),'.mat'));        
                %fname = fullfile(pwd, 'data', 'perturb_ini',...
                %    strcat(fname_str,'-v',int2str(i),'.mat'));
            end
            save(fname)
        end
    end
end
