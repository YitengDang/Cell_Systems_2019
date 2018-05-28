%% Plots normalized noise of p_out averaged over all values of p_in
clear all
close all
clc

%% Parameters
gridsize = 15;
N=gridsize^2;
Con = 16;
K = 8;
a0 = 5;
hill = 2;

%% Load pin-pout data
fname_str = strrep(sprintf('pin_pout_N%d_Con%d_K%d_a0_%d_hill%.2f', N, Con, K, a0*10, hill),...
    '.','p');
fname = fullfile(pwd, 'figures', 'pin_pout', 'data', fname_str);
load(fname, 'prob');

%% 
p = (0:N)/N;
noise = (prob - mean(prob, 1))./std(prob, 0, 1);

figure();
hold on
histogram(noise, 'Normalization', 'pdf');
ksdensity(reshape(noise, 1, numel(noise)));