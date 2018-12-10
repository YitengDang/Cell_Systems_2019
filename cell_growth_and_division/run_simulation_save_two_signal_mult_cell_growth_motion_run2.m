% Simulates the two gene system with cell motility, as modelled by Brownian
% motion of each cell
close all
clear all
% maxNumCompThreads(4);
% warning off

%% (1) input parameters
% lattice parameters
gz = 16;
N = gz^2;
a0 = 1.5;

% movement parameters
sigma_D = 0; %0.1;

% circuit parameters 
M_int = [1 1; -1 0];
Con = [946 870];
Coff = [1 1];
K = [125 803; 261 0];% K(i,j): sensitivity of type i to type j molecules
lambda = [1 1.2]; % diffusion length (normalize first to 1)
lambda12 = lambda(2)/lambda(1);
hill = Inf;
noise = 0;

% growth parameters
sigma_rcell = 0; %0.1; % initial spread in cell radii
rcell = 0.2;
Rcell = rcell*a0;
c_growth = 0; %2; %1.2; % logistic growth rate
%K_growth = 0.15; % carrying capacity (in terms of area covered)
ini_density = 2*pi/sqrt(3)*rcell^2;
k_growth = 1.5; % max. fold-change of density
K_growth = k_growth*ini_density;

% division parameters
r_div_mean = 0.3; % units of 1/a0

% calculate division probability
sigma = 0.1*r_div_mean;
% p_div = cdf('Normal',[0.2 0.25 0.28],r_div_mu,sigma)

% initial conditions
p0 = [0.5 0.5];
iniON = round(p0*N);
I0 = [0.2 0.2];
dI = 0.01;
InitiateI = 1; % 0: no, 1: yes

% generate cell_type (0 case type 1, 1 case type 2)
cell_type = zeros(N,1);

% simulation parameters
tmax = 1000;
mcsteps = 0;
nruns = 35;

save_folder = 'H:\My Documents\Multicellular automaton\temp\Network_19';
%save_folder = 'W:\staff-homes\d\yitengdang\My Documents\Multicellular automaton\temp';
%save_folder = 'W:\staff-homes\d\yitengdang\My Documents\Multicellular automaton\temp';

%{
fname_lbl = strrep(sprintf('N%d_iniON_%d_%d_M_int_%d_%d_%d_%d_a0_%.1f_Con_%d_%d_K_%d_%d_%d_%d_lambda_%.1f_%.1f_mcsteps_%d', ...
    N, iniON(1), iniON(2), M_int(1,1), M_int(1,2), M_int(2,1), M_int(2,2), ...
    a0, Con(1), Con(2), K(1,1), K(1,2), K(2,1), K(2,2),...
    lambda(1), lambda(2), mcsteps), '.', 'p');
%}

% generate initial lattice (dist, pos)
%[dist, pos] = init_dist_hex(gz, gz);

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

% Plot growth rate vs. rcell
%{
h = figure;
density_all = 0.01:0.01:K_growth;
rcell_all = sqrt(density_all/(2*pi/sqrt(3)) );
plot(density_all, c_growth.*density_all.*(1-density_all/K_growth) );
%plot(rcell_all, c_growth.*density_all.*(K_growth-density_all) );
%ylim([0 1]);
%xlabel('$r_{cell}$');
xlabel('Density');
ylabel('logistic growth rate');

%fprintf('max(mu) = %.3f at rho = %.3f \n', K_growth^2/2*(1-1/(2*c_growth)), K_growth/(2*c_growth) );
%}
%% ----------- simulation ------------------------------------
for irun=1:nruns
    
% settings
t = 0;
disp_mol = 12;
showI = 0;
cells_hist = {};
positions_all = {};
rcell_hist = {};
max_rejections_reached = 0;

% Initial position -> initiate on hexagonal lattice
nodisplay = 1;
[pos, dist] = initial_cells_random_markov_periodic(gz, mcsteps, rcell, nodisplay);

% store cell positions
positions_all{end+1} = pos;

% (1) generate initial state
%
%iniON = round(p0*N);
cells = zeros(N, 2);
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
% store cell states
cells_hist{end+1} = cells; %{cells(:, 1), cells(:, 2)};


% Initiate cells on a lattice, with random cell sizes
%rcell_all = zeros(N, 1);
rcell_all = rcell+randn(N, 1)*sigma_rcell*rcell;
rcell_hist{end+1} = rcell_all;
%Rcell_all = rcell_all*a0;

% Plot lattice
%{
hin = figure;
%set(hin, 'Position', [100 100 800 800]);
[h_cells, h_borders, a0_px] = reset_cell_figure_cell_growth(hin, pos, rcell); %reset_cell_figure_cell_growth(hin, pos, rcell_all);
%update_figure_periodic_scatter_cell_growth(h_cells, h_borders, cells, t, disp_mol, showI, a0, dist, rcell_all, a0_px);
update_figure_periodic_cell_motion_cell_sizes(h_cells, h_borders,...
    cells, t, disp_mol, showI, a0, a0_px, dist, pos, rcell_all);
%}
% Update cells
%[cells_out, changed] = update_cells_two_signals_multiply_finite_Hill(cells, dist, M_int, a0,...
%    Rcell, Con, Coff, K, lambda, hill, noise);
[cells_out, changed] = ...
    update_cells_two_signals_multiply_finite_Hill_w_cell_sizes(cells,...
        dist, M_int, a0, a0*rcell_all, Con, Coff, K, lambda, hill, noise);

% update positions
%[pos, dist, rejections] = update_cell_positions(gz, rcell_all, pos, dist, sigma_D);
[pos, dist, rejections] = update_cell_positions_diff_cell_sizes(...
    gz, rcell_all, pos, dist, sigma_D);
%fprintf('Update position rejections = %d \n', rejections);

% update cell sizes
mu_cells = calc_growth_rate(rcell_all, c_growth, K_growth); % calculate density
rcell_all_new = rcell_all*sqrt((1+mu_cells)); % sqrt: growth rate is for biomass

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
    else    
        % Check that distances are within range
        rcell_mat = repmat(rcell_all_new, 1, N) + repmat(rcell_all_new', N, 1); % sum of radii of pairs of cells
        cond_mat = (dist > rcell_mat);
        if all(cond_mat(dist>0))
            % if passed, update cell radii
            rcell_all = rcell_all_new;
            rcell_hist{end+1} = rcell_all;
            % update cell sizes 
            mu_cells = calc_growth_rate(rcell_all, c_growth, K_growth); % calculate growth rate
            rcell_all_new = rcell_all*sqrt((1+mu_cells)); % sqrt: growth rate is for biomass
        else
            % growth hits a maximum, don't update cell radii
            warning('Updating cell radii would lead to overlapping cells!');
            rcell_hist{end+1} = rcell_all;
        end
    end
    
    % Update figure
    %{
    update_figure_periodic_scatter_cell_growth(h_cells, h_borders, cells, t,...
        disp_mol, showI, a0, dist, rcell_all, a0_px);
    update_figure_periodic_cell_motion_cell_sizes(h_cells, h_borders,...
        cells, t, disp_mol, showI, a0, a0_px, dist, pos, rcell_all);
    %}
    
    % Update cell states
    %[cells_out, changed] = update_cells_two_signals_multiply_finite_Hill(cells, dist, M_int, a0,...
    %    Rcell, Con, Coff, K, lambda, hill, noise);
    [cells_out, changed] = ...
        update_cells_two_signals_multiply_finite_Hill_w_cell_sizes(cells,...
            dist, M_int, a0, a0*rcell_all, Con, Coff, K, lambda, hill, noise);
    
    % update positions
    %[pos, dist, rejections] = update_cell_positions(gz, rcell_all, pos, dist, sigma_D);
    [pos, dist, rejections] = update_cell_positions_diff_cell_sizes(...
        gz, rcell_all, pos, dist, sigma_D);
    if rejections >= 10^6
        % some cells are not being updated because the max. number of
        % MC rejections has been reached
        max_rejections_reached = 1;
    end
    %fprintf('Update position rejections = %d \n', rejections);
    
    % continue if cells are still growing and cell states changed
    if round(mu_cells, 4)>0
        cont_sim = 1;
    elseif t_growth_stop==Inf
        t_growth_stop = t;
        fprintf('Growth ceased at t = %d \n', t);
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
ext = '.mat';
%label = '';
% default file name
if InitiateI
    I_ini_str = sprintf('_I_ini_%.2f_%.2f', I0(1), I0(2));
else
    I_ini_str = '';
    I0 = Inf;
end
if max_rejections_reached
    label = '_max_rejections_reached';
else
    label = '';
end
%fname_str = strrep(sprintf('two_signals_growing_cells_t_out_%d',...
%    t_out), '.', 'p');
fname_str = strrep(sprintf('two_signals_rcell_sigma_%.1f_K_growth_%.1f_sigma_D_%.3f_t_out_%d_neg_control',...
    sigma_rcell, k_growth, sigma_D, t_out), '.', 'p');
%fname_str = strrep(sprintf('two_signals_horiz_TW_initial_sigma_rcell_%.1f_K_growth_1p5_ini_density_remains_wave',...
%    sigma_rcell), '.', 'p');
%{
if ~cont_sim
    survival_str = 'SURVIVED';
else
    survival_str = 'destroyed';
end
%}
%fname_str = sprintf('%s_rcell_0p3_sigma_rcell_0_no_growth_%s', ini_state_fname,...
%    survival_str);
%fname_str = sprintf('%s_sigma_rcell_0p1_K_growth_1p5_%s', ini_state_fname,...
%    survival_str);

% check if filename already exists
i=1;
fname = fullfile(save_folder, strcat(fname_str, '-v', num2str(i),...
    label, ext));
while exist(fname, 'file') == 2
    i=i+1;
    fname = fullfile(save_folder, strcat(fname_str, '-v', num2str(i),...
        label, ext));
end

save_vars = {N, a0, K, Con, Coff, M_int, hill, noise,...
    p0, I0, rcell, sigma_rcell, sigma_D, ...
    lambda12, I_ini_str, mcsteps, c_growth, K_growth, r_div_mean};
save_vars_lbl = {'N', 'a0', 'K', 'Con', 'Coff', 'M_int', 'hill', 'noise',...
    'p_ini', 'I_ini', 'rcell', 'sigma_rcell', 'sigma_D',....
    'lambda12','I_ini_str', 'mcsteps', 'c_growth', 'K_growth', 'r_div_mean'};

save_consts_struct = cell2struct(save_vars, save_vars_lbl, 2);
positions = pos;
%positions_all = pos_hist;
distances = dist;
t_out = t;

qsave = 1;
if qsave
    save(fname, 'save_consts_struct', 'cells_hist', 't_out',...
        'changed', 'positions', 'distances', 'positions_all', 'rcell_hist');
    fprintf('Saved simulation: %s ; \n', fname);
end

end
%}
%% Load simulation
%{
%folder = 'H:\My Documents\Multicellular automaton\temp';
folder = 'W:\staff-homes\d\yitengdang\My Documents\Multicellular automaton\temp';
version = 1;
t_out = 247;
fname_str = strrep(sprintf('two_signals_growing_cells_t_out_%d-v%d',...
    t_out, version), '.', 'p');
fname_str = 'two_signals_growing_cells_moving_sigma_D_0p100_t_out_121-v1';

ext = '.mat';
fname = fullfile(folder, strcat(fname_str, ext));
load(fname);
%}
%% Save as movie
%{
[folder, file] = fileparts(fname);
fname_out = fullfile(folder, strcat(file, '.avi'));
frame_rate = 4;

save_movie = 0;
% Options
%frame_rate = 5; % frames/second
format = 'Motion JPEG AVI'; %movie format 
% 'Motion JPEG AVI' <- default, works best
% 'Uncompressed AVI' <- high quality(?), large file
% 'MPEG-4' <- .mp4
% 'Archival' <- unknown ext
% 'Motion JPEG 2000' <- unknown ext

% Initiate movie
if save_movie
    myVideo = VideoWriter(fname_out, format); %, 'Uncompressed AVI');
    myVideo.FrameRate = frame_rate;  % Default 30
    open(myVideo);
end
% replay trajectory externally
%tt = 0;
h = figure;
clf(h, 'reset');
[h_cells, h_borders, a0_px] = reset_cell_figure_cell_growth(h, pos, rcell); %reset_cell_figure_cell_growth(hin, pos, rcell_all);

if isempty(positions_all)
    % same positions every time step
    for tt=1:length(cells_hist)
        pause(0.01);
        cells = cells_hist{tt};
        rcell_all = rcell_hist{tt};
        %update_figure_periodic_scatter_cell_growth(h_cells, h_borders,...
        %    cells, tt, disp_mol, showI, a0, dist, rcell_all, a0_px);
        update_figure_periodic_cell_motion_cell_sizes(h_cells, h_borders,...
            cells, tt, disp_mol, showI, a0, a0_px, dist, pos, rcell_all);
        %update_cell_figure_external(h, pos, cells, cell_type, tt-1, disp_mol, rcell)                
        %frames(t) = getframe(h);
        frame = getframe(h);
        if save_movie
            writeVideo(myVideo, frame);
        end
    end
else
    % new positions every time step
    for tt=1:length(cells_hist)
        pause(0.01);
        cells = cells_hist{tt};
        pos = positions_all{tt};
        rcell_all = rcell_hist{tt};
        %update_figure_periodic_scatter_cell_growth(h_cells, h_borders,...
        %    cells, tt, disp_mol, showI, a0, dist, rcell_all, a0_px);   
        update_figure_periodic_cell_motion_cell_sizes(h_cells, h_borders,...
            cells, tt, disp_mol, showI, a0, a0_px, dist, pos, rcell_all);
        %frames(t) = getframe(h);
        frame = getframe(h);
        if save_movie
            writeVideo(myVideo, frame);
        end
    end
end
%

% add final state to show equilibrium (not applicable if not in
% equilibrium)
%{
cells = cells_hist{t};
t = t+1;
update_cell_figure_external(h, pos, cells, cell_type, t-1, disp_mol, rcell)     
frames(t) = getframe(h);
frame = getframe(h);
writeVideo(myVideo, frame);
%}
if save_movie
    close(myVideo);
end
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
