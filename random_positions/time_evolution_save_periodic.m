%% Time evolution of lattice with visualization
clear all
close all
set(0,'defaulttextinterpreter', 'latex');
%%
% geometric parameters
Lx = 1;
n = 12; % nmax = L/R (square packing)
rcell = 0.2;
ini_alg = 0; % 1: Markov MC, 0: random placement
mcsteps = 10^4; % for Markov MC
%--------------------
N = n^2; % total number of particles
%N = round(eta*L^2/(pi*R^2));
Ly = sqrt(3)/2*Lx;
R = rcell*Lx/n; % disc radius
r_av = sqrt(Lx*Ly/N)/2; % estimate for NND of random distribution (Clark & Evans, 1954)
eta = N*pi*R^2/(Lx*Ly); %packing fraction 
if eta > pi/(2*sqrt(3)) % max. packing fraction
    disp('packing fraction too high! Abort evaluation');
    pause(10);
end

% Circuit Parameters
K_all = [6 12 20 16];
Con_all = [15 15 15 8];
noise = 0;
a0_all = [1.5 1.5 1.5 0.5];
rcell = 0.2;

% K=[6 12 20 16], Con=[15 15 15 8], a0 = [1.5 1.5 1.5 0.5]
%Klist = [3 6 10 15 17 20]; %(a0=1.5) [3 6 10 15 17 20]; % (a0=0.5) [10 10 14 19 16 18]; %[3]; 
%Conlist = [24 21 21 20 14 14]; %(a0=1.5) [24 21 21 20 14 14]; % (a0=0.5) [5 21 16 14 8 6]; %[24];

% initial config
p_all = 0:0.1:1;
ntrials = 10;

for idx_K=1:4
    K = K_all(idx_K);
    Con = Con_all(idx_K);
    a0 = a0_all(idx_K);
    for idx_p=1:numel(p_all)
        p = p_all(idx_p);
        iniON = round(N*p);

        % Save parameters
        save_folder = 'H:\My Documents\Multicellular automaton\data\random_positions';
        labels = {'RandomPlacement', sprintf('MarkovMC_%d_MCsteps', mcsteps)};
        fname_str = strrep(sprintf('N%d_Lx%.2f_a0_%.2f_K%d_Con%d_n%d_%s', N, Lx,...
            a0, K, Con, round(p*N), labels{ini_alg+1}), '.', 'p');

        fprintf('p = %.2f, %d trials to go %d \n', p, ntrials);
        for trial=1:ntrials
            disp(trial);
            cells_hist = {};
            t=0;

            % Initiate lattice
            if ini_alg
                % (1) Markov MC
                [pos, dist, fN0] = initial_cells_random_markov_periodic(n, Lx, R, mcsteps);
            else
                % (2) random placement
                [pos, dist] = initial_cells_random_periodic_alt(n, Lx, Ly, R); % random config
            end

            % display initial lattice
            %hin = figure(1);
            cells = zeros(N, 1);
            cells(randperm(N, iniON)) = 1;
            %update_figure_periodic_long(pos+R, N, Lx, Ly, R, cells, t)
            cells_hist{end+1} = cells;

            % Time evolution
            changed = 1;
            while changed
                % calculate new state
                pause(0.2);
                [cells_out, changed] = ...
                    update_cells_noise(cells, dist_new, Con, K, a0, rcell*a0, noise);

                %pause(1);
                % update cells & figure
                t = t+1;
                cells = cells_out;
                cells_hist{end+1} = cells;
                %update_figure_periodic_long(pos+R, N, Lx, Ly, R, cells, t);
            end

            close all

            % Save trajectory
            v=1;
            fname = fullfile(save_folder, strcat(fname_str, '-v', num2str(v), '.mat') );
            while exist(fname, 'file')==2
                v = v+1;
                fname = fullfile(save_folder, strcat(fname_str, '-v', num2str(v), '.mat') );
            end
            save(fname, 'cells_hist', 'dist', 't', 'N', 'Lx', 'a0', 'mcsteps', 'Con', 'K',...
                'noise', 'rcell', 'p');
            fprintf('Saved file: %s \n', fname);
        end
    end
end