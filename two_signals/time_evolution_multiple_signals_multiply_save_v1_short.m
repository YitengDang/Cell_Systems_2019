close all
clear all
maxNumCompThreads(4);
%warning off
%% (1) input parameters
%
% lattice parameters
gz = 15;
N = gz^2;
a0 = 1.5;
rcell = 0.2;
Rcell = rcell*a0;

% circuit parameters
Con = [18 16];
Coff = [1 1];
M_int = [0 1; -1 1];
K = [0 25; 11 4]; % K(i,j): sensitivity of type i to type j molecules
lambda = [1 1.2]; % diffusion length (normalize first to 1)
hill = Inf;
noise = 0;

% initial conditions
%p0 = [0.2 0.6];
%iniON = round(p0*N);
I0 = [0 0];
dI = 0.01;
InitiateI = 0; % 0: no, 1: yes

% generate cell_type (0 case type 1, 1 case type 2)
cell_type = zeros(N,1);

% simulation parameters
tmax = 10000;
nruns = 80;

% pos, dist
[dist, pos] = init_dist_hex(gz, gz);
%[pos,ex,ey] = init_cellpos_hex(gridsize,gridsize);
%dist = dist_mat(pos,gridsize,gridsize,ex,ey);

%{
fname_str = strrep(sprintf('N%d_iniON_%d_%d_M_int_%d_%d_%d_%d_a0_%.1f_Con_%d_%d_K_%d_%d_%d_%d_lambda_%.1f_%.1f', ...
    N, iniON(1), iniON(2), M_int(1,1), M_int(1,2), M_int(2,1), M_int(2,2), ...
    a0, Con(1), Con(2), K(1,1), K(1,2), K(2,1), K(2,2),...
    lambda(1), lambda(2)), '.', 'p');
%}
% check parameters
idx = (M_int == 0);
if ~all(K(idx)==0)
    fprintf('K has wrong entries! \n');
    warning('K has wrong entries!');
end

%}
%{
% TO DO: vectorize
dist_vec = a0*dist(1,:);
r = dist_vec(dist_vec>0); % exclude self influence
fN1 = sum(sinh(Rcell)*sum(exp((Rcell-r)./lambda(1)).*(lambda(1)./r)) ); % calculate signaling strength
fN2 = sum(sinh(Rcell)*sum(exp((Rcell-r)./lambda(2)).*(lambda(2)./r)) ); % calculate signaling strength

% nearest neighbour interaction strength
fprintf('activator fij(a0) = %.4f \n', sinh(Rcell)*sum(exp((Rcell-a0)./lambda(1)).*(lambda(1)./a0)))
fprintf('inhibitor fij(a0) = %.4f \n', sinh(Rcell)*sum(exp((Rcell-a0)./lambda(2)).*(lambda(2)./a0)))
%}

%% (2) Load parameters from saved trajectory
%{
% with parameters saved as structure array 
% load data
%data_folder = 'H:\My Documents\Multicellular automaton\app\Multicellularity-2.1\data\time_evolution';
data_folder = 'D:\Multicellularity\app\Multicellularity-2.1\data\time_evolution';
[file, path] = uigetfile(fullfile(data_folder, '\*.mat'), 'Load saved simulation');
load(fullfile(path, file));

s = save_consts_struct;
N = s.N;
a0 = s.a0;
K = s.K;
Con = s.Con;
Coff = s.Coff;
M_int = s.M_int;
hill = s.hill;
noise = s.noise;
rcell = s.rcell;
cells = cells_hist{1};
lambda = [1 s.lambda12];
%p0 = s.p_ini;
tmax =  s.tmax;
gz = sqrt(N);
Rcell = rcell*a0;
[dist, pos] = init_dist_hex(gz, gz);

% simulation parameters
tmax = 10^4;
nruns = 10;
cell_type = zeros(N,1);
 
% Initial I
InitiateI = 0;
I0 = [0 0];
s_fields = fieldnames(s);
for i=1:numel(s_fields)
    if strcmp(s_fields{i},'I_ini_str')
        if ~isempty(s.I_ini_str)
            I0 = s.I_ini;
            InitiateI = 1;
        end
    end
end

%}
% Check whether loaded trajectory is same as simulation
%{
eq = 2*ones(numel(cells_hist_2), 1);
for i=1:numel(cells_hist_2)
    eq(i) = all(all(cells_hist{i} == cells_hist_2{i}));
end
%}

%% Loop over p_ini
% -- setup
[p1, p2] = meshgrid(0:0.1:1, 0:0.1:1);
%p1 = s.p_ini(1);
%p2 = s.p_ini(2);
% folder to save in
%folder = fullfile('H:\My Documents\Multicellular automaton\data\two_signals\time_evolution',...
%    'temp');
folder = fullfile('D:\Multicellularity\data\two_signals\time_evolution',...
    'vs_pini_batch3');

% default file name
sim_ID = 'two_signal_mult';
M_int_str = sprintf('M_int%d_%d_%d_%d', M_int(1,1), M_int(1,2), M_int(2,1), M_int(2,2));
K_str = sprintf('K%.1f_%.1f_%.1f_%.1f', K(1,1), K(1,2), K(2,1), K(2,2));
Con_str = sprintf('Con%.1f_%.1f', Con(1), Con(2));
Coff_str = sprintf('Coff%.1f_%.1f', Coff(1), Coff(2));
I_ini_str = '';
if InitiateI
    I_ini_str = sprintf('_I1_I2_%.2f_%.2f', I0(1), I0(2));
end
lambda12 = lambda(2);
    
%% -- loop
for i1 = 1:size(p1, 1)
    for i2 = 1:size(p1, 2)
        fprintf('p1 = %.1f, p2 = %.1f \n', p1(i1, i2),  p2(i1, i2) )
        p0 = [ p1(i1, i2)  p2(i1, i2) ];
        
        % default file name
        p_ini_str = sprintf('p_ini%.2f_%.2f', p0(1), p0(2));
        fname_str = strrep(sprintf('%s_%s_hill%.2f_N%d_a0_%.2f_%s_%s_%s_noise%.1f_%s%s_l12_%.1f_rcell%.1f_tmax%d',...
            sim_ID, M_int_str, hill, N, a0, K_str, Con_str, Coff_str, noise,...
            p_ini_str, I_ini_str, lambda12, rcell, tmax), '.', 'p');
        
        % skip parameters that are complete
        if exist(fullfile(folder, sprintf('%s-v%d.mat', fname_str, nruns)), 'file')==2
           disp('parameter set done, continued')
           continue
        end        
        
        % Dynamics
        for run=1:nruns
            % skip runs that are complete
            if exist(fullfile(folder, sprintf('%s-v%d.mat', fname_str, run)), 'file')==2
               fprintf('run %d exists, continued \n', run);
               continue
            end
            
            disp(run);
            cells_hist = {};
            %t_out = 0;
            %changed = 1;

            % generate initial lattice
            iniON = round(p0*N);
            cells = zeros(N, 2);
            for i=1:numel(iniON)
                cells(randperm(N,iniON(i)), i) = 1;
                if InitiateI && hill==Inf
                    fprintf('Generating lattice with I%d(t=0)... \n', i);
                    [cells_temp, test, I_ini] = generate_I_new(cells(:, i), I0(i), I0(i)+dI, dist, a0);
                    cells(:,i) = cells_temp;
                    fprintf('Generated initial I%d: %.2f; in range (Y=1/N=0)? %d; \n', i, I_ini, test);
                end
            end

            % store initial config
            cells_hist{end+1} = cells; %{cells(:, 1), cells(:, 2)};
            
            %-------------dynamics-----------------------------------------
            t = 0;
            period = Inf; %default values
            t_onset = Inf; 
            [cellsOut, changed] = update_cells_two_signals_multiply_finite_Hill(cells, dist, M_int, a0,...
                    Rcell, Con, Coff, K, lambda, hill, noise);
            
            while changed && t < tmax && period==Inf
                %pause(0.2);
                t = t+1;
                cells = cellsOut;
                cells_hist{end+1} = cells; %{cells(:, 1), cells(:, 2)};
                [period, t_onset] = periodicity_test_short(cells_hist); 

                %update_cell_figure_continuum(app, pos, dist, a0, cells, app.Time, cell_type, disp_mol, 0);
                [cellsOut, changed] = update_cells_two_signals_multiply_finite_Hill(cells, dist, M_int, a0,...
                    Rcell, Con, Coff, K, lambda, hill, noise);
            end
            t_out = t; % save final time
            
            if period==Inf
                fprintf('t_out = %d, no periodicity found \n', t_out);
            end
            %--------------------------------------------------------------
            % Save result
            if isempty(cells_hist)
                warning('No simulation data to save!');
            end

            ext = '.mat';
            label = '';

            if exist(folder, 'dir') ~= 7
                warning('Folder does not exist! ');
            end

            % filename 
            %filename = strrep(sprintf('%s_%s_hill%.2f_N%d_a0_%.2f_%s_%s_%s_noise%.1f_%s%s_l12_%.1f_rcell%.1f_tmax%d',...
            %    sim_ID, M_int_str, hill, N, a0, K_str, Con_str, Coff_str, noise, p_ini_str, I_ini_str, lambda12, rcell, tmax), '.', 'p');

            % check if filename already exists
            i=1;
            fname = fullfile(folder, strcat(fname_str, '-v', num2str(i), label, ext));
            while exist(fname, 'file') == 2
                i=i+1;
                fname = fullfile(folder, strcat(fname_str, '-v', num2str(i), label, ext));
            end

            save_vars = {N, a0, K, Con, Coff, M_int, hill, noise, p0, I0, rcell,...
                tmax, lambda12, sim_ID, I_ini_str};
            save_vars_lbl = {'N', 'a0', 'K', 'Con', 'Coff', 'M_int', 'hill', 'noise', 'p_ini', 'I_ini', 'rcell',...
                'tmax', 'lambda12', 'sim_ID', 'I_ini_str'};
            save_consts_struct = cell2struct(save_vars, save_vars_lbl, 2);

            save(fname, 'save_consts_struct', 'cells_hist', 't_out',...
                'changed', 'period', 't_onset');
            fprintf('Saved simulation: %s ; \n', fname);

        end
%--------------end for loop for run----------------------------------------
    end
end