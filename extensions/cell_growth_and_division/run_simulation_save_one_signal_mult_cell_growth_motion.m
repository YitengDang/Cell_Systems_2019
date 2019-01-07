% Simulates the two gene system with cell motility, as modelled by Brownian
% motion of each cell
close all
clear all
% maxNumCompThreads(4);
% warning off

%% (1) input parameters
% lattice parameters
gz = 15;
N = gz^2;
a0 = 0.5;

% movement parameters
sigma_D = 0.1;

% circuit parameters 
M_int = 1;
Con = 8;
Coff = 1;
K = 15;
hill = Inf;
noise = 0;

% growth parameters
sigma_rcell = 0; %0.1; % initial variation in cell sizes
rcell = 0.2;
Rcell = rcell*a0;
c_growth = 0; %2; % max. growth rate
ini_density = 2*pi/sqrt(3)*rcell^2;
k_growth = 1.5; % max. fold-change of density
K_growth = k_growth*ini_density; % carrying capacity (in terms of area covered)

% division parameters
r_div_mean = 0.3; % units of 1/a0

% calculate division probability
sigma = 0.1*r_div_mean;
% p_div = cdf('Normal',[0.2 0.25 0.28],r_div_mu,sigma)

% initial conditions
p0_all = (0:N)/N;

%p0 = 0.5;
%iniON = round(p0*N);

InitiateI = 0; % 0: no, 1: yes
I0 = 0.4;
dI = 0.01;

% generate cell_type (0 case type 1, 1 case type 2)
cell_type = zeros(N,1);

% simulation parameters
tmax = 100;
mcsteps = 10^4;
%prec = 8;
nruns = 20;
display_sim = 1;

% save_folder = 'H:\My Documents\Multicellular automaton\temp';
% save_folder = 'H:\My Documents\Multicellular automaton\temp';
%save_folder = 'W:\staff-bulk\tnw\BN\HY\Shared\Yiteng\two_signals\moving_cells_one_signal\random_position_non_moving';
%save_folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\moving_cells_one_signal\random_position_non_moving';
save_folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\moving_cells_one_signal\random_position_D_0p1';
% save_folder = 'D:\temp';

%% ----------- simulation ------------------------------------
for irun=1:nruns
    for idx_p0=1:numel(p0_all)
        p0 = p0_all(idx_p0);
        iniON = round(p0*N);
    
        fprintf('Run %d \n', irun);

        % settings
        t = 0;
        disp_mol = 1;
        showI = 0;
        cells_hist = {};
        positions_all = {};
        rcell_hist = {};
        %rcell_all = [];
        
        % (1) generate initial state
        %
        %iniON = round(p0*N);
        cells = zeros(N, 1);
        for i=1:numel(iniON)
            cells(randperm(N,iniON(i)), i) = 1;
            if InitiateI && hill==Inf
                fprintf('Generating lattice with I%d(t=0)... \n', i);
                dI = 0.01;
                [cells_temp, test, I_ini] = generate_I_new(cells(:, i), I0(i), I0(i)+dI, dist, a0);
                cells(:,i) = cells_temp;
                fprintf('Generated initial I%d: %.2f; in range (Y=1/N=0)? %d; \n', i, I_ini, test);
            end
        end
        %}
        % (2) load initial state from input
        %{
        signal_count = 2;
        folder = 'D:\Multicellularity\app\data\system_states';
        [status, cells, ini_state_fname] = manual_input_state(signal_count, folder, N);
        %}

        % Initiate cells with random cell sizes
        %rcell_all = zeros(N, 1);
        rcell_all = rcell+randn(N, 1)*sigma_rcell*rcell;
        %Rcell_all = rcell_all*a0;

        % generate initial lattice (dist, pos) (! do after initiating cell radii)
        %[dist, pos] = init_dist_hex(gz, gz);
        nodisplay = 0;
        [pos, dist] = initial_cells_random_markov_periodic_diff_cell_sizes(gz, mcsteps, rcell_all, nodisplay);

        % store data 
        cells_hist{end+1} = cells; % cell states {cells(:, 1), cells(:, 2)};
        positions_all{end+1} = pos; %  cell positions
        rcell_hist{end+1} = rcell_all;
        %%
        %------------
        % Plot lattice
        hin = figure;
        %set(hin, 'Position', [100 100 800 800]);
        [h_cells, h_borders, a0_px] = reset_cell_figure_cell_growth(hin, pos, rcell); %reset_cell_figure_cell_growth(hin, pos, rcell_all);
        %update_figure_periodic_scatter_cell_growth(h_cells, h_borders, cells, t, disp_mol, showI, a0, dist, rcell_all, a0_px);
        update_figure_periodic_cell_motion_cell_sizes(h_cells, h_borders,...
            cells, t, disp_mol, showI, a0, a0_px, dist, pos, rcell_all);
        %%
        %------------
        % Update cells
        [cells_out, changed] = ...
            update_cells_noise_hill_diff_cell_sizes(cells, dist, M_int, Con, K, a0,...
            rcell_all*a0, noise, hill);

        % update positions
        %[pos, dist, rejections] = update_cell_positions(gz, rcell_all, pos, dist, sigma_D);
        [pos, dist, rejections] = update_cell_positions_diff_cell_sizes(...
            gz, rcell_all, pos, dist, sigma_D);
        %fprintf('Update position rejections = %d \n', rejections);

        % update cell sizes
        %mu_cells = calc_growth_rate(rcell_all, c_growth, K_growth); % calculate density
        %rcell_all_new = rcell_all*sqrt((1+mu_cells)); % sqrt: growth rate is for biomass
        [rcell_all_new, delta_A_cells] = calc_cell_sizes_logistic_growth(rcell_all,...
            c_growth, K_growth);
        %%
        %------------
        cont_sim = 1;
        period_ub = Inf; % upper bound for periodicity (from periodicity_test_short)
        t_growth_stop = Inf; % time at which growth ceases (because of saturation)
        while t<tmax && cont_sim
            disp(t);
            t = t+1;
            pause(0.1);

            % Store cell states
            cells = cells_out;
            cells_hist{end+1} = cells; %{cells(:, 1), cells(:, 2)};

            % store cell positions
            positions_all{end+1} = pos;

            if t_growth_stop<Inf
                % growth has ceased, don't update cell radii
                rcell_hist{end+1} = rcell_all;
                delta_A_cells = 0;
            else    
                % Check that distances are within range
                rcell_mat = repmat(rcell_all_new, 1, N) + repmat(rcell_all_new', N, 1); % sum of radii of pairs of cells
                cond_mat = (dist > rcell_mat);
                if all(cond_mat(dist>0))
                    % if passed, update cell radii
                    rcell_all = rcell_all_new;
                    rcell_hist{end+1} = rcell_all;
                    % update cell sizes 
                    %mu_cells = calc_growth_rate(rcell_all, c_growth, K_growth); % calculate growth rate
                    %rcell_all_new = rcell_all*sqrt((1+mu_cells)); % sqrt: growth rate is for biomass
                    [rcell_all_new, delta_A_cells] = calc_cell_sizes_logistic_growth(rcell_all,...
                        c_growth, K_growth);
                else
                    % growth hits a maximum, don't update cell radii
                    warning('Updating cell radii would lead to overlapping cells!');
                    rcell_hist{end+1} = rcell_all;
                    delta_A_cells = 0;
                end
            end

            % Update figure
            %update_figure_periodic_scatter_cell_growth(h_cells, h_borders, cells, t,...
            %    disp_mol, showI, a0, dist, rcell_all, a0_px);
            update_figure_periodic_cell_motion_cell_sizes(h_cells, h_borders,...
                cells, t, disp_mol, showI, a0, a0_px, dist, pos, rcell_all);

            % Update cell states
            % Update cells
            [cells_out, changed] = ...
                update_cells_noise_hill_diff_cell_sizes(cells, dist, M_int, Con, K, a0,...
                rcell_all*a0, noise, hill);

            % update positions
            %[pos, dist, rejections] = update_cell_positions(gz, rcell_all, pos, dist, sigma_D);
            [pos, dist, rejections] = update_cell_positions_diff_cell_sizes(...
                gz, rcell_all, pos, dist, sigma_D);
            fprintf('Update position rejections = %d \n', rejections);

            % continue if cells are still growing and cell states changed
            %fprintf('Summed change in area: %.4f \n', sum(abs(delta_A_cells)) );
            if delta_A_cells > 10^(-4)
                cont_sim = 1;
            elseif t_growth_stop==Inf
                t_growth_stop = t;
                fprintf('Growth ceased at t = %d \n', t);
            elseif sigma_D>0
                % always continue if cells are moving
                cont_sim = 1;
            else
                cont_sim = changed;
                % check for periodicity
                [period_ub, ~] = periodicity_test_short(cells_hist(t_growth_stop+1:end)); 
                if period_ub<Inf % stop when periodicity has been found
                    cont_sim = 0;
                end
            end
            % cont_sim = ~round(mu_cells, 3)==0;
        end
        t_out = t;

        % if periodicity found, refine check to find period
        if period_ub<Inf
            t_check_init = t_growth_stop;
            decimals = Inf;
            [period, t_onset] = periodicity_test_detailed(cells_hist, t_check_init,...
                period_ub, decimals);
            t_out = t_onset + period; 
        end
        %% Save trajectory
        % default file name
        if InitiateI
            I_ini_str = sprintf('_I_ini_%.2f', I0);
        else
            I_ini_str = '';
            I0 = Inf;
        end

        fname_str = strrep(sprintf('one_signal_M_int_%d_mcsteps_%d_sigma_D_%.3f_no_growth_tmax_%d_iniON_%d_t_out_%d',...
            M_int, mcsteps, sigma_D, tmax, iniON, t_out), '.', 'p');
        ext = '.mat';
        label = '';

        % check if filename already exists
        i=1;
        fname = fullfile(save_folder, strcat(fname_str, '-v', num2str(i), label, ext));
        while exist(fname, 'file') == 2
            i=i+1;
            fname = fullfile(save_folder, strcat(fname_str, '-v', num2str(i), label, ext));
        end

        save_vars = {N, a0, K, Con, Coff, M_int, hill, noise, p0, I0, rcell,...
            sigma_rcell, I_ini_str, mcsteps, c_growth, K_growth, r_div_mean};
        save_vars_lbl = {'N', 'a0', 'K', 'Con', 'Coff', 'M_int', 'hill', 'noise',...
            'p_ini', 'I_ini', 'rcell', 'sigma_rcell', 'I_ini_str', 'mcsteps',...
            'c_growth', 'K_growth', 'r_div_mean'};

        save_consts_struct = cell2struct(save_vars, save_vars_lbl, 2);
        positions = pos;
        distances = dist;

        qsave = 1;
        if qsave
            save(fname, 'save_consts_struct', 'cells_hist', 't_out',...
                'changed', 'positions_all', 'positions', 'distances', 'rcell_hist');
                fprintf('Saved simulation: %s ; \n', fname);
        end
    end
end

%% Plot cell mass density vs. rcell
%{
h = figure;
rcell_all = 0.01:0.01:0.5;
density_all = 2*pi/sqrt(3)*rcell_all.^2;
plot(rcell_all, density_all);
ylim([0 1]);
xlabel('$r_{cell}$');
ylabel('fraction area covered');
%}
%{
% Plot growth rate vs. rcell
h = figure;
density_all = 0.01:0.01:K_growth;
rcell_all = sqrt(density_all/(2*pi/sqrt(3)) );
plot(rcell_all, mu_max.*density_all.*(K_growth-density_all) );
%ylim([0 1]);
xlabel('$r_{cell}$');
ylabel('logistic growth rate');
%}

%% Simulate growth only
%{
% Initiate cells on a lattice, with random cell sizes
%rcell_all = zeros(N, 1);
rcell_all = rcell+randn(N, 1)*0.1*rcell;
cells = zeros(N, 2);

% generate initial lattice
%[dist, pos] = init_dist_hex(gz, gz);
nodisplay = 1;
mcsteps = 0;
[pos, dist] = initial_cells_random_markov_periodic(gz, mcsteps, rcell, nodisplay);

% Plot lattice
hin = figure;
set(hin, 'Position', [100 100 800 800]);
[h_cells, h_borders, a0_px] = reset_cell_figure_cell_growth(hin, pos, rcell); %reset_cell_figure_cell_growth(hin, pos, rcell_all);
plot_handle = h_cells;
pause(0.1);
t = 0;
disp_mol = 12;
showI = 0;
update_figure_periodic_scatter_cell_growth(h_cells, h_borders, cells, t, disp_mol, showI, a0, dist, rcell_all, a0_px);

% calculate density
A_system = sqrt(3)/2; % area of system
A_cells = pi*mean(rcell_all.^2);
rho = A_cells/A_system;
mu_cells = mu_max*rho*(K_growth-rho);

% update cell sizes
rcell_all_new = rcell_all*sqrt((1+mu_cells)); % sqrt: growth rate is for biomass
    
tmax = 100;
cont_sim = 1;
while t<tmax && cont_sim
    disp(t);
    t = t+1;
    pause(1);
    
    % Check that distances are within range
    rcell_mat = repmat(rcell_all_new, 1, N) + repmat(rcell_all_new', N, 1); % sum of radii of pairs of cells
    cond_mat = (dist > rcell_mat);
    if ~all(cond_mat(dist>0))
        % update cell radii so they don't overlap anymore
        warning('overlapping!');
        cont_sim = 0;
        break
    end
    % if passed, update cell radii
    rcell_all = rcell_all_new;
    
    % Update figure
    update_figure_periodic_scatter_cell_growth(h_cells, h_borders, cells, t, disp_mol, showI, a0, dist, rcell_all, a0_px);
    
    % calculate density
    A_system = sqrt(3)/2; % area of system
    A_cells = pi*mean(rcell_all.^2);
    rho = A_cells/A_system;
    mu_cells = mu_max*rho*(K_growth-rho); % logistic growth, no noise
    disp(mu_cells);
    
    % update cell sizes
    rcell_all_new = rcell_all*sqrt((1+mu_cells)); % sqrt: growth rate is for biomass
   
    % stop if cells don't grow anymore
    cont_sim = ~round(mu_cells, 3)==0;
end
%}

%% (2) Load parameters from saved trajectory
%{
% with parameters saved as structure array 
% load data
data_folder = 'H:\My Documents\Multicellular automaton\app\Multicellularity-2.1\data\time_evolution';
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
%tmax = 100;
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