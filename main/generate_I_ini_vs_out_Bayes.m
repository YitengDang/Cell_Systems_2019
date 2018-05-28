% This script uses the calculated P(p_out, I_out | p_in, I_in) to perform
% Bayesian inference and calculate P(p_in, I_in | p_out, I_out) maps
close all
clear variables
%warning off

% Parameters of the system
gridsize = 15;
N = gridsize^2;
a0 = 0.5;
Rcell = 0.2*a0;
K = 15;
Con = 8;
%a0all = [0.5 1.5 1.5];
%Kall = [15 6 20];
%Conall = [8 15 15];

p1 = floor(0.1*N)/N;
dp = round(0.05*N)/N;
p_all = p1:dp:(1-p1);
I_all = -0.05:0.05:0.4;

% use hexagonal lattice
[dist, pos] = init_dist_hex(gridsize, gridsize);

% filename of the saved data
fname_str = strrep(sprintf('pIin_pIout_N%d_Con_%d_K_%d_a0_%.1f', ...
    N, Con, K, a0), '.', 'p');

% Calculate the signaling strength
dist_vec = dist(1,:);
r = a0*dist_vec(dist_vec>0); % exclude self influence
fN = sum(sinh(Rcell)*sum(exp(Rcell-r)./r)); % calculate signaling strength

%% load the calculated map
%uigetfile(fullfile(pwd, 'data', 'pIin_pIout'));
v = 1;
fname = fullfile(pwd, 'data', 'pIin_pIout', strcat(fname_str,'-v', int2str(v), '.mat'));
load(fname)