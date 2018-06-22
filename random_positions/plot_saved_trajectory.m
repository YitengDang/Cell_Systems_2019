%% Plots surviving cells vs initial density
clear all
close all
clc
%%
% geometric parameters
L = 2;
R = 0.02; % disc radius

% circuit parameters
K = 16;
S = 8;

% data properties 
n = 9; % nmax = L/R (square packing)
N = n^2; % total number of particles

r_av = sqrt(L^2/N)/2; % estimate for NND of random distribution (Clark & Evans, 1954)
eta = N*pi*R^2/L^2; %packing fraction 

straux = '(\d+)';
fpattern = sprintf('N%d_K%d_S%d_L%.2f_R%.2f_n%d-v%s', N, K, S, L, R, n, straux); 

[pos, dist, rej] = initial_cells_random(N, L, R);

%% Load data

path = fullfile(pwd, 'data');

% Get all file names in the directory
listing = dir(path);
num_files = numel(listing)-2; %first two entries are not useful
count = 0;
for i = 1:num_files
    filename = listing(i+2).name;
    % remove extension and do not include txt files
    [~,name,ext] = fileparts(filename);
    if strcmp(ext, '.mat')
        count = count + 1;
        names{count} = name;
    end
end

% Store fraction of final cells
idx=0; %counts number of loaded simulations
for i = 1:numel(names)
    % first get the filename given by the student
    [tokens, ~] = regexp(names{i},fpattern,'tokens','match');
    if numel(tokens) > 0
        disp(names{i}) % displays the file name
        idx = idx+1;
        % load the data
        load(fullfile(path,strcat(names{i},'.mat')), 'cells_hist', 't');
        hin=figure(1);
        for i=0:t
            update_figure(hin, pos, N, L, R, cells_hist{i+1});
        end
        k=waitforbuttonpress;
    end
end