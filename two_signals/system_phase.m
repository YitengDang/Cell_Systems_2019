%% determines the phase the system is in for each of the component molecules
clear all
close all
%% Parameters

% (1) input parameters
% lattice parameters
gz = 15;
N = gz^2;
a0 = 1.5;
rcell = 0.2;
Rcell = rcell*a0;

% circuit parameters
Con = [18 16];
Coff = [1 1];
M_int = [1 1; -1 -1];
K = [3 12; 13 20]; % K(i,j): sensitivity of type i to type j molecules
lambda = [1 1.2]; % diffusion length (normalize first to 1)
hill = Inf;
noise = 0;

% pos, dist
[dist, pos] = init_dist_hex(gz, gz);

% check parameters
idx = (M_int == 0);
if ~all(K(idx)==0)
    fprintf('K has wrong entries! \n');
    warning('K has wrong entries!');
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
%[dist, pos] = init_dist_hex(gz, gz);
%}

%% Calculate fN
% TO DO: vectorize
dist_vec = a0*dist(1,:);
r = dist_vec(dist_vec>0); % exclude self influence
fN1 = sum(sinh(Rcell)*sum(exp((Rcell-r)./lambda(1)).*(lambda(1)./r)) ); % calculate signaling strength
fN2 = sum(sinh(Rcell)*sum(exp((Rcell-r)./lambda(2)).*(lambda(2)./r)) ); % calculate signaling strength

% nearest neighbour interaction strength
fprintf('activator fij(a0) = %.4f \n', sinh(Rcell)*sum(exp((Rcell-a0)./lambda(1)).*(lambda(1)./a0)))
fprintf('inhibitor fij(a0) = %.4f \n', sinh(Rcell)*sum(exp((Rcell-a0)./lambda(2)).*(lambda(2)./a0)))
%}


%% Calculate phase
num_mol = 2; % number of molecules
cond = zeros(num_mol, 4); % test 4 conditions per molecule
fN = [fN1 fN2];

all_ON_all = zeros(num_mol,num_mol);
all_OFF_all = zeros(num_mol,num_mol);
ON_remains_ON_all = zeros(num_mol,num_mol);
OFF_remains_OFF_all = zeros(num_mol,num_mol);

% determine gij for each interaction
for i=1:num_mol
    for j=1:num_mol
        all_ON_all(i,j) = gij(1, i, j, M_int, fN, K, Coff, Con, 1);
        all_OFF_all(i,j) = gij(0, i, j, M_int, fN, K, Coff, Con, 1);
        ON_remains_ON_all(i,j) = gij(1, i, j, M_int, fN, K, Coff, Con, -1);
        OFF_remains_OFF_all(i,j) = gij(0, i, j, M_int, fN, K, Coff, Con, -1);
    end
end

% final result for each gene
all_ON = prod(all_ON_all, 2);
all_OFF = prod(all_OFF_all, 2);
ON_remains_ON = prod(ON_remains_ON_all, 2);
OFF_remains_OFF = prod(OFF_remains_OFF_all, 2);
%%
condition = {'All ON', 'All OFF', 'ON -> ON or OFF -> ON', 'OFF remains OFF'};
molecule_number = (1:num_mol)';
t = table(molecule_number, all_ON, all_OFF, ON_remains_ON, OFF_remains_OFF)

% gij(0, 2, 2, M_int, fN, K, Coff, Con, -1)
function gij_out = gij(gij_target, i, j, M_int, fN, K, Coff, Con, cases)
    % cases = 1 for all ON/OFF, cases = -1 for ON cannot turn OFF / reverse
    if M_int(i,j)==0 % no interaction
        gij_out = 1;
        return
    end
    C_all = [Coff(j); Con(j)];
    
    Gamma = M_int(i,j)*(2*gij_target-1);
    beta = -Gamma; 
    idx = (beta==-1) + 2*(beta==1);
    Cb = C_all(idx);
    alpha = cases*beta;
    idx2 = (alpha==-1) + 2*(alpha==1);
    Ca = C_all(idx2);
    %fprintf('Gamma = %.1f \n', Gamma);
    %fprintf('Ca = %.1f \n', Ca);
    %fprintf('Cb = %.1f \n', Cb);
    %fprintf('beta = %d, idx = %d \n', beta, idx);
    gij_out = (Gamma*(Ca + fN(j)*Cb - K(i,j)) > 0);
end
%{
% Old version: does not work

function out = gij_zero_strong(i, j, M_int, fN, K, Con, Coff, onoff)
    % onoff: (0) all ON (1) all OFF
    S_nei = (1-onoff)*Coff(j)+ onoff*Con(j);
    Y_self = S_nei; % (1-onoff)*Coff(j)+ onoff*Con(j);
    if M_int(i,j)==0
        out = true;
    elseif M_int(i,j)==1
        out = Y_self + fN(j)*S_nei - K(i,j) > 0;
    elseif M_int(i,j)==-1
        out = Y_self + fN(j)*S_nei - K(i,j) < 0;
    end
end

function out = gij_zero_weak(i, j, M_int, fN, K, Con, Coff, onoff)
    % out=0 means test failed, out=1 means success
    % onoff: test whether (0) ON cannot turn OFF (1) OFF cannot turn ON
    % weak_cond: for ON cannot turn OFF and vice versa cases
    S_nei = (1-onoff)*Coff(j)+ onoff*Con(j);
    Y_self = (1-onoff)*Con(j)+ onoff*Coff(j);
    if M_int(i,j)==0
        out = true;
    elseif M_int(i,j)==1
        out = Y_self + fN(j)*S_nei - K(i,j) > 0;
    elseif M_int(i,j)==-1
        out = Y_self + fN(j)*S_nei - K(i,j) < 0;
    end
end
%}
%% Plot phase diagrams
