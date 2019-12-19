% Time evolution of a system without noise and without visualization
close all
clear all
warning off
%%
% Lattice parameters
gz = 15;
N = gz^2;
%a0 = 5.5;
%Rcell = 0.2*a0;
rcell = 0.2;

% loop parameters
a0_list = [0.5 1.5];
Con_list = 1:30;
K_list = 1:40;

% circuit parameters
hill = 2; % Hill coefficient
noise = 0;

% simulation parameters
num_sims = 10;
qsave = 1;
initialID = 'uniform';
sim_ID = 'one_signal_finite_hill';
tmax = 10^3;
prec = 8;

% generate initial lattice
%[dist, pos] = init_dist_hex(gz, gz);
mcsteps = 0;
nodisplay = 1;
[pos, dist] = initial_cells_random_markov_periodic(gz, mcsteps, rcell, nodisplay);
InitiateI = 0;
p_ini = Inf;
I_ini = Inf;
cells_ini = [];
%% Calculate # simulations to do
sim_to_do = zeros(numel(a0_list), numel(Con_list), numel(K_list));
for idx_1=1:numel(a0_list)
    for idx_2=1:numel(Con_list)
        for idx_3=1:numel(K_list)
            this_a0 = a0_list(idx_1);
            this_Con = Con_list(idx_2);
            this_K = K_list(idx_3);

            % Count how many simulations have already been done
            parent_folder = fullfile('N:\tnw\BN\HY\Shared\Yiteng\one_signal_finite_Hill\dynamics',...
                '2019-04-25_vs_K_Con_a0_0p50_1p50_hill2p00_uniform');
            subfolder = strrep(sprintf('a0_%.2f', this_a0), '.', 'p');
            subfolder2 = sprintf('Con_%d_K_%d', this_Con, this_K);
            save_folder = fullfile(parent_folder, subfolder, subfolder2);

            % subfolder
            folder = save_folder;
            if exist(folder, 'dir') ~= 7
                mkdir(folder);
            end
            %}
            if exist(folder, 'dir') ~= 7
                warning('Folder does not exist! ');
                break
            end

            % Filename pattern 
            %--------------------------------------------------------------
            pattern = strrep(sprintf('N%d_a0_%.2f_K_%d_Con_%d_hill_%.2f_t_out_%s-v%s',...
                    N, this_a0, this_K, this_Con, hill,...
                    '(\d+)', '(\d+)'), '.', 'p');
            %--------------------------------------------------------------
            
            listing = dir(folder);
            num_files = numel(listing)-2;
            names = {};
            filecount = 0;
            for i = 1:num_files
                filename = listing(i+2).name;
                % remove extension and do not include txt files
                [~,name,ext] = fileparts(filename);
                if strcmp(ext, '.mat')
                    match = regexp(name, pattern, 'match');
                    %disp(match);
                    if ~isempty(match)
                        filecount = filecount + 1;
                        names{end+1} = name;
                    end
                end
            end

            fprintf('a0 = %.2f, Con = %d, K = %d, sim to do: %d \n', this_a0, this_Con, this_K, num_sims-filecount);
            sim_to_do(idx_1, idx_2, idx_3) = num_sims-filecount;
        end
    end
end
%% Run simulations
for idx_a0=1:numel(a0_list)
    a0 = a0_list(idx_a0);
    for i=1:numel(Con_list)
        Con = Con_list(i);
        for j=1:numel(K_list)
            K = K_list(j);
            fprintf('a0 = %.2f, K=%d, Con=%d \n', a0, K, Con);
            
            subfolder = strrep(sprintf('a0_%.2f', a0), '.', 'p');
            subfolder2 = sprintf('Con_%d_K_%d', Con, K);
            save_folder = fullfile(parent_folder, subfolder, subfolder2);
            
            % Set parameters
            Rcell = rcell*a0;
            dist_vec = a0*dist(1,:);
            r = dist_vec(dist_vec>0); % exclude self influence
            fN = sum(sinh(Rcell)*sum(exp(Rcell-r)./r)); % calculate signaling strength

            % generate cell_type (0 case type 1, 1 case type 2)
            cell_type = zeros(N,1); % all the same here
            
            % default filename
            fname_str_template = strrep(sprintf('N%d_a0_%.2f_K_%d_Con_%d_hill_%.2f',...
            	N, a0, K, Con, hill), '.', 'p');
            display_fig = 0;

            for tests = 1:sim_to_do(idx_a0, i, j)
                time_evolution_save_func_one_signal(...
                    N, a0, Rcell, hill, noise, K, Con, prec, ...
                    dist, pos, sim_ID, mcsteps, InitiateI, p_ini, I_ini, cells_ini, tmax,...
                    save_folder, fname_str_template, display_fig);
            end
        end
    end
end

%{
%% Plot <Xi>(t)
figure();
hold on
plot(0:t, xmean, 'bo-');
plot(0:t, xstd, 'ro-');

%% Plot energy
figure();
plot(0:t, mom); %h = H/N
%}
