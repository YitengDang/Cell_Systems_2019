%{
clear all
close all
clc
%% lattice parameters
gz = 5;
a0 = 1.5;
rcell = 0.2;

% circuit parameters 
M_int = [0 -1; 1 1];
Con = [18 16];
Coff = [1 1];
K = [0 25; 5 10];% K(i,j): sensitivity of type i to type j molecules
lambda = [1 1.2]; % diffusion length (normalize first to 1)
hill = Inf;
noise = 0;

% looping parameters
K12_all = [5 35];
K21_all = [5 35];
K22_all = [5 10 15];
%%
entropy_two_signals_sim_make_dir(gz, a0, rcell, M_int, Con, K, lambda, hill, noise, K12_all, K21_all, K22_all)
%}
function entropy_two_signals_sim_make_dir(parent_folder, gz, a0, rcell, M_int, Con, K, lambda, hill, noise)
    % Constructs set of subfolders for the simulations of the multicellular
    % entropy
    N = gz^2;
    
    M_int_s = sprintf('%d_%d_%d_%d', M_int(1,1), M_int(1,2), M_int(2,1), M_int(2,2));
    a0_s = sprintf('%.2f', a0);
    R_s = sprintf('%.2f', rcell);
    %K_s = sprintf('%d_%d_%d_%d', K(1,1), K(1,2), K(2,1), K(2,2));
    Con_s = sprintf('%d_%d', Con(1), Con(2));
    lambda_s = sprintf('%.1f', lambda(2));
    
    fname_str_combined = strrep(sprintf(...
        'N%d_M_int_%s_a0_%s_rcell_%s_lambda12_%s_Con_%s_hill%.1f_noise%.1f',...
         N, M_int_s, a0_s, R_s, lambda_s, Con_s, hill, noise), '.', 'p');
    %% Make parent folder
    
    if exist(parent_folder, 'dir')~=7
        mkdir(fullfile(parent_folder, fname_str_combined));
    end
    %% Make K folder
    K_s = sprintf('K_%d_%d_%d_%d', ...
        K(1,1), K(1,2), K(2,1), K(2,2));
    folder = fullfile(parent_folder, fname_str_combined, K_s);
    if exist(folder, 'dir')~=7
        mkdir(folder);
    end
end