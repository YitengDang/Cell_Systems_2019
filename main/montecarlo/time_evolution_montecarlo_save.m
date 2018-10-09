%% Plot trajectories (simulated, Langevin) on top of a map of Peq in p, I space
clear all
close all
clc
set(0, 'defaulttextinterpreter', 'latex');

%% Parameters
% Loop parameters
K_all = 7:9;
Con_all = 12:14;

% Fixed parameters
gridsize = 15;
N = gridsize^2;
a0 = 1.5;
rcell = 0.2;
Rcell = rcell*a0;
K = 7;
Con = 12;
alpha = 0; % noise

% load fN, gN
[dist, pos] = init_dist_hex(gridsize, gridsize);
dist_vec = dist(1,:);
r = a0*dist_vec(dist_vec>0); % exclude self influence
fN = sum(sinh(Rcell)*sum(exp(Rcell-r)./r)); % calculate signaling strengthN
gN = sum(sum((sinh(Rcell)*exp(Rcell-r)./r).^2)); % calculate signaling strength

% Simulation parameters
p0 = 0.4;
I0 = 0;
n_runs = 10;
tmax = 1000;
noSpatialOrder = 1;

% for saving
if noSpatialOrder
    I_s = 'no_Theta';
else
    I_s = 'with_Theta';
end
fname_str = strrep(sprintf(...
    'MC_EOM_%s_repression_N%d_a0_%.2f_K_%d_Con_%d_noise_%.1f_p_ini_%.1f_I_ini_%.1f_t_%d',...
    I_s, N, a0, K, Con, alpha, p0, I0, tmax), '.', 'p');
parent_folder = 'N:\tnw\BN\HY\Shared\Yiteng\one_signal_repression';
%%
for idxK = 1:numel(K_all)
    for idxCon = 1:numel(Con_all)
        subfolder = strrep(sprintf('K_%.1f_Con_%.1f', K, Con), '.', 'p');
        folder = fullfile(parent_folder, subfolder);
        if exist(folder, 'dir') ~= 7
            mkdir(folder);
            fprintf('Made directory', folder);
        end

        % Single simulation 
        p_t = [];
        %I_t = [];
        Theta_t = [];

        p = p0;
        Theta = fN*((2*p-1)^2 + 4*p*(1-p)*I0);

        t = 0;
        p_t(end+1) = p;
        Theta_t(end+1) = Theta;
        %[theta, p, pe] = update_montecarlo(theta, p, N, Con, K, fN, gN, alpha);
        [Theta, p, pe] = update_montecarlo_repression(Theta, p, N, Con, K, fN, gN, alpha, noSpatialOrder);
        while rand > pe && t < tmax
            t = t+1;
            p_t(end+1) = p;
            Theta_t(end+1) = Theta;
            %[theta, p, pe] = update_montecarlo(theta, p, N, Con, K, fN, gN, alpha);
            [Theta, p, pe] = update_montecarlo_repression(Theta, p, N, Con, K, fN, gN, alpha, noSpatialOrder);
        end

        % save trajectory
        qsave = 1;
        if qsave
            i = 1;
            fname = fullfile(folder, strcat(fname_str,'-v',int2str(i), '.mat'));
            while exist(fname, 'file') == 2
                i=i+1;
                fname = fullfile(folder, ...
                    strcat(fname_str,'-v',int2str(i),'.mat'));
            end
            save(fname);
        end
    end
end