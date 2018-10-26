% Predicts whether a travelling wave can propagate according to nearest
% neighbour interactions
clear all
close all
set(0,'defaulttextinterpreter', 'latex')
%% Parameters
% Manual input
gz = 15;
N = gz^2;
a0 = 1.5;
rcell = 0.2;
Rcell = rcell*a0;
lambda = [1 1.2];
M_int = [1 1; -1 0]; % network 19 
%{
M_int = [0 1; -1 1]; % network 15 reversed
M_int = [1 -1; 1 0]; % network 15
M_int = [1 1; -1 0]; % network 19
M_int = [1 -1; 1 1]; % network 33
M_int = [-1 -1; 1 1]; % network 34
M_int = [-1 1; -1 1]; % network 36
%}
Con = [18 16];
K = zeros(2);
%K = [0 9; 11 4];

% get pos, dist
mcsteps = 0;
[pos, dist] = initial_cells_random_markov_periodic(gz, mcsteps, rcell);

% Obtain from simulation
folder = 'L:\BN\HY\Shared\Yiteng\two_signals\batch_sim_all_topologies_run2\selected';
subfolder = 'patterns\Network 33';
fname_str = 'tripple_wave_diagonally';
fname = fullfile(folder, subfolder, fname_str);
load(fname, 'save_consts_struct');

Con = save_consts_struct.Con;
K = save_consts_struct.K;
N = save_consts_struct.N;
gz = sqrt(N);
a0 = save_consts_struct.a0;
rcell = save_consts_struct.rcell;
Rcell = rcell*a0;
lambda = save_consts_struct.lambda;
M_int = save_consts_struct.M_int;

%%
% calculate fN
idx = gz + round(gz/2); % pick cell not at corner of grid
dist_vec = a0*dist(idx,:);
r = dist_vec(dist_vec>0);
fN = zeros(2,1);
fN(1) = sum(sinh(Rcell)*exp((Rcell-r)./lambda(1)).*(lambda(1)./r));
fN(2) = sum(sinh(Rcell)*exp((Rcell-r)./lambda(2)).*(lambda(2)./r));

% calculate fnn
fnn = zeros(1, 2);
rnn = a0;
fnn(1) = sinh(Rcell)*exp((Rcell-rnn)./lambda(1)).*(lambda(1)./rnn) ; % calculate signaling strength
fnn(2) = sinh(Rcell)*exp((Rcell-rnn)./lambda(2)).*(lambda(2)./rnn) ; % calculate signaling strength

fnnn = zeros(2,2);
rnnn = [sqrt(3) 2].*a0;
fnnn(1,:) = sinh(Rcell)*exp((Rcell-rnnn)./lambda(1)).*(lambda(1)./rnnn);
fnnn(2,:) = sinh(Rcell)*exp((Rcell-rnnn)./lambda(2)).*(lambda(2)./rnnn);

% save folder 
save_folder = 'H:\My Documents\Multicellular automaton\figures\trav_wave_stability';
fname_str_default = strrep(sprintf('N%d_a0_%.1f_rcell_%.1f_lambda12_%.1f_Con_%d_%d_K_%d_%d_%d_%d',...
    N, a0, rcell, lambda(2), Con(1), Con(2), K(1,1), K(1,2), K(2,1), K(2,2)), '.', 'p');

%% Calculate Y_all
% specify wave type and characteristics
wave_types_str = {'straight', 'inward bend', 'outward bend'};
wave_type = 1;
num_waves = 2; % number of waves
bandwidth = 1; % width of band of cells of one type

% order of states: F, M, B, E
states_perm = [4 2 1 3]; 
%{
states_perm = [2 4 3 1]; % network 15
states_perm = [4 3 1 2]; % network 19
states_perm = [3 4 2 1]; % network 33
states_perm = [4 2 1 3]; % network 33/34
states_perm = [2 4 3 1]; % network 36
%}
default_states = [0 0; 0 1; 1 0; 1 1];

% calculate Y(alpha)
% neighbour data
switch wave_type
    case 1
        n_nei = [2	0	0	4; % EF
            2	2	0	2; % F
            2	2	2	0; % M 
            0	2	2	2; % B
            0	0	2	4; % EB
            0	0	0	6]; % E
    case 2
        n_nei = [3	0	0	3;
            2	3	0	1;
            1	2	3	0;
            0	1	2	3;
            0	0	1	5;
            0	0	0	6];
    case 3
        n_nei = [1	0	0	5;
            2	1	0	3;
            3	2	1	0;
            0	3	2	1;
            0	0	3	3;
            0	0	0	6];
end
states = default_states(states_perm, :); % F, M, B, E
types = [states(4,:); states(1,:); states(2,:); states(3,:); states(4,:); states(4,:)];
Y_self = (Con-1).*types + 1;

% calculate Y_nei
Y_nei = zeros(size(n_nei, 1), 2);
for i=1:6
    Y_nei(i,1) = fnn(1)*n_nei(i,:)*((Con(1)-1)*states(:,1)+1);
    Y_nei(i,2) = fnn(2)*n_nei(i,:)*((Con(2)-1)*states(:,2)+1);
end
%
% estimate p
tmp = num_waves*bandwidth;
tmp2 = [tmp tmp tmp gz-3*tmp];
p = (tmp2*states)/gz;

% calculate Y_mf (mean-field contribution)
z = 6; % coordination number
Y_mf = zeros(1, 2);
Y_mf(1) = (fN(1) - z*fnn(1))*( Con(1)*p(1) + (1-p(1)) );
Y_mf(2) = (fN(2) - z*fnn(2))*( Con(2)*p(2) + (1-p(2)) );
%disp(Y_mf)

Y_all = Y_self + Y_nei + Y_mf;

%% Check conditions for specific parameter set
conditions_met = zeros(6,1);
% (1) E_F?F
% (2) F?M
% (3) M?B
% (4) B?E
% (5) E_B?E
% (6) E?E 
targets = [1 2 3 4 4 4]; % recall F, M, B, E
output_state_list = zeros(6, 2);
for i=1:6
    %i = 1;
    Y = Y_all(i, :); % 2x1
    output_state = prod(((Y - K).*M_int>0) + (1-abs(M_int)), 2)';
    target_state = states(targets(i),:);
    conditions_met(i) = all(output_state==target_state);
    output_state_list(i,:) = output_state;
end
%disp(output_state_list);
disp(conditions_met);

%% Check for a set of parameters
var_x_all = 1:10:1000; %variable to plot on x axis
var_y_all = 1:10:1000; %variable to plot on y axis
var_z_all = 1:10:1000; %variable to plot on z axis
%K_22_all = 1:10;
%trav_wave_cond_met = zeros(numel(K_12_all), numel(K_21_all), numel(K_22_all));
trav_wave_cond_met = zeros(numel(var_x_all), numel(var_y_all), numel(var_z_all));

for k=1:numel(var_z_all)
    for j=1:numel(var_y_all)
        for i=1:numel(var_x_all)
            thisK = K;
            thisK(1,2) = var_x_all(i); %K12
            thisK(2,1) = var_y_all(j); %K21
            %thisK(2,2) = K_22_all(k);
            thisK(1,1) = var_z_all(k); %K11
            
            conditions_met = zeros(6,1);
            % (1) E_F?F
            % (2) F?M
            % (3) M?B
            % (4) B?E
            % (5) E_B?E
            % (6) E?E 
            targets = [1 2 3 4 4 4]; % recall F, M, B, E
            output_state_list = zeros(6, 2);
            for i2=1:6
                %i = 1;
                Y = Y_all(i2, :); % 2x1
                output_state = prod(((Y - thisK).*M_int>0) + (1-abs(M_int)), 2)';
                target_state = states(targets(i2),:);
                conditions_met(i2) = all(output_state==target_state);
                output_state_list(i2,:) = output_state;
            end
            %disp(output_state_list);
            %disp(conditions_met);
            
            trav_wave_cond_met(i,j,k) = all(conditions_met);
        end        
    end
end

%% 3D scatter plot
%[K_21_mesh, K_12_mesh, K_22_mesh] = meshgrid(K_21_all, K_12_all, K_22_all);
[K_y_mesh, K_x_mesh, K_z_mesh] = meshgrid(var_y_all, var_x_all, var_z_all); %N.B. x,y order changed by meshgrid

h = figure;
hold on
grid on
colors = {'b', 'r', 'g'};
idx1 = squeeze(trav_wave_cond_met==1);
%s1 = scatter3(K_12_mesh(idx1), K_21_mesh(idx1), K_22_mesh(idx1),...
%    150, 'o', 'filled', 'MarkerFaceColor', colors{wave_type});
s1 = scatter3(K_x_mesh(idx1), K_y_mesh(idx1), K_z_mesh(idx1),...
    150, 'o', 'filled', 'MarkerFaceColor', colors{wave_type});
s1.MarkerFaceAlpha = 0.5;
xlim([var_x_all(1)-1 var_x_all(end)+1]);
ylim([var_y_all(1)-1 var_y_all(end)+1]);
zlim([var_z_all(1)-1 var_z_all(end)+1]);
az = -45; el = 20;
view(az, el);
xlabel('$$K^{(12)}$$');
ylabel('$$K^{(21)}$$');
zlabel('$$K^{(22)}$$');
set(gca, 'FontSize', 20);
set(h, 'Units', 'inches', 'position', [1 1 10 6]);
daspect([1 1 1])

qsave = 0;
fname_str = strrep(sprintf('Trav_wave_stability_analytical_v2_scatter3_%s_num_waves_%d',...
    wave_types_str{wave_type}, num_waves), ' ', '_');    
fname = fullfile(save_folder, fname_str);
save_figure(h, 10, 6, fname, '.pdf', qsave);
%}

%% Save data
%{
save_folder = 'H:\My Documents\Multicellular automaton\figures\trav_wave_stability\data';
fname_str_data = sprintf('Trav_wave_stability_analytical_v2_%s_num_waves_%d',...
    fname_str_default, num_waves);  
fname = fullfile(save_folder, fname_str_data);
save(fname, 'K_12_all', 'K_21_all', 'K_22_all', 'trav_wave_cond_met', ...
    'num_waves', 'bandwidth', 'wave_type');
%}

%% Load data 
%{
trav_wave_cond_met_2 = trav_wave_cond_met;
%
load_folder = 'H:\My Documents\Multicellular automaton\figures\trav_wave_stability\data';
%fname_str_data = sprintf('Trav_wave_stability_analytical_%s', fname_str_default); 
fname_str_data = 'Trav_wave_stability_analytical_N225_a0_1p5_rcell_0p2_lambda12_1p2_Con_18_16_K_0_9_11_1';
fname = fullfile(load_folder, fname_str_data);
load(fname, 'K_12_all', 'K_21_all', 'K_22_all', 'trav_wave_cond_met');

compare = (squeeze(trav_wave_cond_met(wave_type,:,:,:)) == trav_wave_cond_met_2);
all(compare(:))
%}