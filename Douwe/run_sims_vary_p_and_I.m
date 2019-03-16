%% Setting scope simulation
maxNumCompThreads(6); % limits the number of cores used by MATLAB

tmax = 10000;
% Folder for storing simulations
%parent_folder = 'D:\Documenten\Studie\Master\Internship Delft\test';
parent_folder = 'K:\bn\hy\Shared\Douwe\Recent Results\Vary_p_and_I';

% Specify number of simulations

% p range:
p1_all = 0.1;%0:0.1:1;
p2_all = 0.1:0.2:0.5;%0:0.1:1;

% I range
I1_all = 0;%-0.1:0.1:0.7;
I2_all = -0.1:0.1:0.7;

% Number of simulations for each unique (p1,p2,I1,I2):
n_sims = 39;

% Genetic Circuit for which simulations are done
network = 33;

% Base simulation.mat file name
base_name = sprintf('Network_%g',network);

% Calc number of simulations
total_sims = length(p1_all)*length(p2_all)*length(I1_all)*length(I2_all)*n_sims;
%% Loading and setting parameters
% Simulation file containing parameters should be only file param subfolder 
% of parent_folder
param_folder = sprintf('%s%s',parent_folder,'/param');
cd(param_folder)
file_name = dir('*.mat');
file_name = file_name.name;
load(file_name)

% Load parameters
N = save_consts_struct.N;
a0 = save_consts_struct.a0;
lambda12 = save_consts_struct.lambda12;
rcell = save_consts_struct.rcell;
M_int = save_consts_struct.M_int; 
Con = save_consts_struct.Con;
K = save_consts_struct.K;
lambda = [1 lambda12];
Coff = [1 1];
Rcell = rcell*a0;
gz = sqrt(N);

hill = Inf;
noise = 0;

InitiateI = 1;
cell_type = zeros(N,1);

mcsteps = 0;
growth_rate = 0;
R_division = 0;
sigma_D = 0;
%% Run Simulations
counter = 0;
for i = 1:n_sims
    
    for p1 = p1_all
        save_folderp1 = fullfile(parent_folder, strrep(sprintf('\\p1_%g',p1),'.','p'));
        
        if ~exist(save_folderp1,'dir')
            mkdir(save_folderp1)
        end
        
        for p2 = p2_all
            save_folderp2 = fullfile(save_folderp1, strrep(sprintf('\\p2_%g',p2),'.','p'));
            
            if ~exist(save_folderp2,'dir')
                mkdir(save_folderp2)
            end
            
            
            for I1 = I1_all
                save_folder_I1 = fullfile(save_folderp2, strrep(sprintf('\\I1_%.1f',I1),'.','p'));
                
                if ~exist(save_folder_I1,'dir')
                    mkdir(save_folder_I1)
                end
                
                
                for I2 = I2_all
                    p0 = [p1 p2];
                    I0 = [I1 I2];
                    sim_ID = base_name;
                    
                    save_folder = fullfile(save_folder_I1, strrep(sprintf('\\I2_%.1f',I2),'.','p'));
                    
                    if ~exist(save_folder,'dir')
                        mkdir(save_folder)
                    end
               
                    [cells_hist, period, t_onset] = time_evolution_save_func_efficient_checks_I_p(...
                        N, a0, Rcell, lambda, hill, noise, M_int, K, Con, Coff,...
                        distances, positions, sim_ID, mcsteps, InitiateI, p0, I0, tmax, save_folder,i);
                    
                    %print info to screen
                    counter = counter + 1;
                    fprintf('%g percent completed\n',counter/total_sims * 100)
                end
            end
        end
    end
end