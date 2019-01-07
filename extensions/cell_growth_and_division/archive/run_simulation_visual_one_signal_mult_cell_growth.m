% Simulates the two gene system with cell motility, as modelled by Brownian
% motion of each cell
% WARNING: INCORRECT logistic growth function 

close all
clear all
% maxNumCompThreads(4);
% warning off

%% (1) input parameters
% lattice parameters
gz = 18;
N = gz^2;
a0 = 1.5;

% movement parameters
sigma_D = 0;

% circuit parameters 
M_int = 1;
Con = 16;
Coff = 1;
K = 18;
hill = Inf;
noise = 0;

% growth parameters
rcell = 0.15;
Rcell = rcell*a0;
mu_max = 1.2; % max. growth rate
K_growth = 0.7; % carrying capacity (in terms of area covered)

% division parameters
r_div_mean = 0.3; % units of 1/a0

% calculate division probability
sigma = 0.1*r_div_mean;
% p_div = cdf('Normal',[0.2 0.25 0.28],r_div_mu,sigma)

% initial conditions
p0 = 0.5;
iniON = round(p0*N);
I0 = 0.4;
dI = 0.01;
InitiateI = 0; % 0: no, 1: yes

% generate cell_type (0 case type 1, 1 case type 2)
cell_type = zeros(N,1);

% simulation parameters
tmax = 50;
mcsteps = 0;
prec = 8;
% settings
disp_mol = 1;
showI = 0;
cells_hist = {};
pos_hist = {};
rcell_hist = {};

%{
fname_lbl = strrep(sprintf('N%d_iniON_%d_%d_M_int_%d_%d_%d_%d_a0_%.1f_Con_%d_%d_K_%d_%d_%d_%d_lambda_%.1f_%.1f_mcsteps_%d', ...
    N, iniON(1), iniON(2), M_int(1,1), M_int(1,2), M_int(2,1), M_int(2,2), ...
    a0, Con(1), Con(2), K(1,1), K(1,2), K(2,1), K(2,2),...
    lambda(1), lambda(2), mcsteps), '.', 'p');
%}

% generate initial lattice (dist, pos)
%[dist, pos] = init_dist_hex(gz, gz);
nodisplay = 1;
[pos, dist] = initial_cells_random_markov_periodic(gz, mcsteps, rcell, nodisplay);

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
%
h = figure;
rcell_all = 0.01:0.01:0.5;
density_all = 2*pi/sqrt(3)*rcell_all.^2;
plot(rcell_all, density_all);
ylim([0 1]);
xlabel('$r_{cell}$');
ylabel('fraction area covered');
%}

% Plot growth rate vs. rcell
h = figure;
density_all = 0.01:0.01:K_growth;
rcell_all = sqrt(density_all/(2*pi/sqrt(3)) );
plot(rcell_all, mu_max.*density_all.*(K_growth-density_all) );
%ylim([0 1]);
xlabel('$r_{cell}$');
ylabel('logistic growth rate');

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

%% ----------- simulation ------------------------------------
t = 0;
% generate initial state
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
% store cell states
cells_hist{end+1} = cells; %{cells(:, 1), cells(:, 2)};

% Initiate cells on a lattice, with random cell sizes
%rcell_all = zeros(N, 1);
rcell_all = rcell+randn(N, 1)*0.1*rcell;
rcell_hist{end+1} = rcell_all;
%Rcell_all = rcell_all * a0;

% Plot lattice
hin = figure;
set(hin, 'Position', [100 100 800 800]);
[h_cells, h_borders, a0_px] = reset_cell_figure_cell_growth(hin, pos, rcell); %reset_cell_figure_cell_growth(hin, pos, rcell_all);
update_figure_periodic_scatter_cell_growth(h_cells, h_borders, cells, t, disp_mol, showI, a0, dist, rcell_all, a0_px);

% Update cells
[cells_out, changed] = ...
    update_cells_noise_hill_diff_cell_sizes(cells, dist, M_int, Con, K, a0,...
    rcell_all*a0, noise, hill, prec);

% update cell sizes
mu_cells = calc_growth_rate(rcell_all, mu_max, K_growth); % calculate density
rcell_all_new = rcell_all*sqrt((1+mu_cells)); % sqrt: growth rate is for biomass

cont_sim = 1;
while t<tmax && cont_sim
    disp(t);
    t = t+1;
    pause(0.2);
    
    % Store cell states
    cells = cells_out;
    cells_hist{end+1} = cells; %{cells(:, 1), cells(:, 2)};

    % Check that distances are within range
    rcell_mat = repmat(rcell_all_new, 1, N) + repmat(rcell_all_new', N, 1); % sum of radii of pairs of cells
    cond_mat = (dist > rcell_mat);
    if all(cond_mat(dist>0))
        % if passed, update cell radii
        rcell_all = rcell_all_new;
        rcell_hist{end+1} = rcell_all;
        % update cell sizes 
        mu_cells = calc_growth_rate(rcell_all, mu_max, K_growth); % calculate growth rate
        rcell_all_new = rcell_all*sqrt((1+mu_cells)); % sqrt: growth rate is for biomass
    else
        % don't update cell radii
        warning('overlapping!');
        rcell_hist{end+1} = rcell_all;
    end
    
    % Update figure
    update_figure_periodic_scatter_cell_growth(h_cells, h_borders, cells, t,...
        disp_mol, showI, a0, dist, rcell_all, a0_px);
    
    % Update cell states
    [cells_out, changed] = ...
        update_cells_noise_hill_diff_cell_sizes(cells, dist, M_int, Con, K, a0,...
        rcell_all*a0, noise, hill, prec);

    % stop if cells don't grow anymore
    % cont_sim = ~round(mu_cells, 3)==0;
end
t_out = t;

%% Save trajectory
save_folder = 'H:\My Documents\Multicellular automaton\temp';

% default file name
if InitiateI
    I_ini_str = sprintf('_I_ini_%.2f', I0);
else
    I_ini_str = '';
    I0 = Inf;
end

fname_str = strrep(sprintf('one_signal_growing_cells_M_int_%d_t_out_%d',...
    M_int, t_out), '.', 'p');
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
    I_ini_str, mcsteps, mu_max, K_growth, r_div_mean};
save_vars_lbl = {'N', 'a0', 'K', 'Con', 'Coff', 'M_int', 'hill', 'noise', 'p_ini', 'I_ini', 'rcell',...
    'I_ini_str', 'mcsteps', 'mu_max', 'K_growth', 'r_div_mean'};

save_consts_struct = cell2struct(save_vars, save_vars_lbl, 2);
positions = pos;
distances = dist;

qsave = 1;
if qsave
    save(fname, 'save_consts_struct', 'cells_hist', 't_out',...
        'changed', 'positions', 'distances', 'rcell_hist');
        fprintf('Saved simulation: %s ; \n', fname);
end
%}
%% Load simulation
%{
folder = 'H:\My Documents\Multicellular automaton\temp';
version = 1;
t_out = 247;
fname_str = strrep(sprintf('two_signals_growing_cells_t_out_%d-v%d',...
    t_out, version), '.', 'p');
%fname_str = 'one_signal_growing_cells_t_out_247_trav_wave-v1';
ext = '.mat';
fname = fullfile(folder, strcat(fname_str, ext));
load(fname);
%}
%% Save as movie
[folder, file] = fileparts(fname);
fname_out = fullfile(folder, strcat(file, '.avi'));
frame_rate = 5;

save_movie = 1;
if save_movie
    % Options
    %frame_rate = 5; % frames/second
    format = 'Motion JPEG AVI'; %movie format 
    % 'Motion JPEG AVI' <- default, works best
    % 'Uncompressed AVI' <- high quality(?), large file
    % 'MPEG-4' <- .mp4
    % 'Archival' <- unknown ext
    % 'Motion JPEG 2000' <- unknown ext

    % Initiate movie
    myVideo = VideoWriter(fname_out, format); %, 'Uncompressed AVI');
    myVideo.FrameRate = frame_rate;  % Default 30
    open(myVideo);

    % replay trajectory externally
    %tt = 0;
    h = figure;
    clf(h, 'reset');
    [h_cells, h_borders, a0_px] = reset_cell_figure_cell_growth(h, pos, rcell); %reset_cell_figure_cell_growth(hin, pos, rcell_all);
    if isempty(pos_hist)
        % same positions every time step
        for tt=1:length(cells_hist)
            cells = cells_hist{tt};
            rcell_all = rcell_hist{tt};
            update_figure_periodic_scatter_cell_growth(h_cells, h_borders,...
                cells, tt-1, disp_mol, showI, a0, dist, rcell_all, a0_px);
            %update_cell_figure_external(h, pos, cells, cell_type, tt-1, disp_mol, rcell)                
            %frames(t) = getframe(h);
            frame = getframe(h);
            writeVideo(myVideo, frame);
        end
    else
        % new positions every time step
        for tt=1:length(cells_hist)
            cells = cells_hist{tt};
            pos = pos_hist{tt};
            rcell_all = rcell_hist{tt};
            update_figure_periodic_scatter_cell_growth(h_cells, h_borders,...
                cells, tt-1, disp_mol, showI, a0, dist, rcell_all, a0_px);               
            %frames(t) = getframe(h);
            frame = getframe(h);
            writeVideo(myVideo, frame);
        end
    end

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

    close(myVideo);

end