%% Plots a map of Peq in p, I space
% Calculates Peq from contributions of both genes separately
clear all
close all
clc
warning off
%% Parameters
% lattice parameters
gz = 10;
N = gz^2;
a0 = 1.5;
rcell = 0.2;
Rcell = rcell*a0;

% circuit parameters 
M_int = [0 -1; 1 1];
Con = [18 16];
Coff = [1 1];
K = [0 25; 5 10];% K(i,j): sensitivity of type i to type j molecules
lambda = [1 1.2]; % diffusion length (normalize first to 1)
hill = Inf;
noise = 0;

% initial conditions
p0 = [0.5 0.5];
iniON = round(p0*N);
I0 = [0 0];
dI = 0.01;
InitiateI = 0; % 0: no, 1: yes

% generate cell_type (0 case type 1, 1 case type 2)
cell_type = zeros(N,1);

% load fN, gN (for perfect lattice)
[dist, pos] = init_dist_hex(gz, gz);
dist_vec = dist(1,:);
r = a0*dist_vec(dist_vec>0); % exclude self influence
fN = sum(sinh(Rcell)*sum(exp(Rcell-r)./r)); % calculate signaling strengthN
gN = sum(sum((sinh(Rcell)*exp(Rcell-r)./r).^2)); % calculate signaling strength

% Alternative arrangement
% mcsteps = 0;
% [pos, dist] = initial_cells_random_markov_periodic(gz, mcsteps, rcell);

%% Calculate P_eq
% Version where we consider the two genes separately
p = (0:N)/N;
%I = 0;
I = linspace(0, 1, 100);
%[pm, Im] = meshgrid(p,I);

Peq = zeros(numel(p), numel(p), numel(I), numel(I));
ip1 = 1; 
iI1 = 1;
%for ip1=1:N+1
    for ip2=1:N+1
        %for iI1=1:numel(I)
            for iI2=1:numel(I)
                fprintf('%d %d %d %d \n', ip1, ip2, iI1, iI2);
                %this_n = [ip1-1 ip2-1];
                %this_I = [I(iI1) I(iI2)];
                this_n = [N ip2-1];
                this_I = [0 iI2];
                
                %[~, ~, pEn] = transition_prob_two_signals_withI(...
                %    this_n, this_I, N, M_int, a0, fN, gN, rcell, K, Con, Coff);
                pEn = Peq_two_signals_withI(this_n, this_I, N, M_int, a0, fN, gN, rcell, K, Con, Coff);
                
                Peq(ip1, ip2, iI1, iI2) = pEn;
            end
        %end
    end
%end

%% Plot result
% Fix I1, I2, plot against p1, p2
%iI1 = 1;
%iI2 = 1;

h=figure(1);
%imagesc(I, p, Peq(:, :, iI1, iI2) );
imagesc(p, I, squeeze(Peq(1, :, 1, :))' );
set(gca, 'YDir', 'Normal');
c = colorbar;
%xlabel('p');
%ylabel('I');
xlabel('$$p^{(2)}$$')
ylabel('$$p^{(1)}$$')
ylabel(c, 'P_{eq}');
caxis([0 1]);
set(gca, 'FontSize', 24);
set(gcf, 'Units','Inches', 'Position', [0 0 10 8]);

%% Equation of motion
%{
%% Calculate delta, sigma (of EOM)
[delta, sigma] = EOM_calc_delta_sigma_two_signals(N, K, Con, fN, gN);

%% Plot gradient vector field
% Multiplicative noise in drift of I
%
%{
% old approach (only N)
Ns = 6^2; %optimal N
fmax = 5; %max amplification
a = Ns*log(fmax);
f = fmax.*exp(-a./N);
f = 1;
%}

% Load fitted surface phi(N, a0)
%fname = 'H:\My Documents\Multicellular automaton\rebuttal\I_correction\I_correction_phi_N_a0_K6_Con15.mat';
%load(fname, 'f');
%phi = f(gridsize, a0);
% phi = maxI Langevin / maxI simulation
% In update rule below, enter 1/phi instead of phi!

% define gradient
delh_delp = @(p,I) (Con-1)*2*fN*(I-1).*(2*p - 1) -(Con+1)*(fN+1)+2*K;
delh_delI = @(p,I) 2*fN*(Con - 1)*p.*(p - 1);

figure(h);
hold on
pv = (0:n:N)/N;
Iv = -0.15:0.08:1;
[pm2, Im2] = meshgrid(pv, Iv);
vf_x = -delh_delp(pm2, Im2);
vf_y = -delh_delI(pm2, Im2);
quiver(pm2, Im2, vf_x, vf_y, 'LineWidth', 1, 'AutoScaleFactor', 1.2,...
    'Color', [0.5 0.5 0.5]);
%}
%% Load simulated trajectories (cellular automaton)
%----Settings--------------------------------------------------------------
% p values of trajectories to plot
p_ini = 0:0.1:1; %(0:5:N)/N;

% load all files from a folder
path = 'H:\My Documents\Multicellular automaton\data\random_positions';
straux = '(\d+)';
labels = {'RandomPlacement', sprintf('MarkovMC_%d_MCsteps', mcsteps)};
choice = 1; %choice of label (randomization procedure)  
%--------------------------------------------------------------------------
listing = dir(path);
num_files = numel(listing)-2; %first two entries are not useful
% variables 
I_ini = zeros(numel(p_ini), 1);
p_out = zeros(numel(p_ini), 1); % store p_out values
eq_times1 = zeros(numel(p_ini), 1);

nruns = 0;
for idx_p = 1:numel(p_ini)
    p = p_ini(idx_p);
    
    % filename pattern (regexp notation)
    fpattern = strrep(sprintf('N%d_Lx%.2f_a0_%.2f_K%d_Con%d_n%d_%s-v%s', N, Lx,...
            a0, K, Con, round(p*N), labels{choice}, straux), '.', 'p');
    %fpattern = strrep(sprintf('N%d_Lx%.2f_n%d_%s-v%s', N, Lx,...
    %     round(p*N), labels{1}, straux), '.', 'p');

    for i = 1:num_files
        filename = listing(i+2).name;
        % remove extension and do not include txt files
        [~,name,ext] = fileparts(filename);
        if strcmp(ext, '.mat')
            [tokens, ~] = regexp(name, fpattern, 'tokens', 'match');
            if numel(tokens) > 0
                disp(name);
                nruns = nruns + 1;
                %names{nruns} = name;
                load(fullfile(path, name), 'cells_hist', 'dist', 'a0_new', 't');

                Non = zeros(t+1, 1);
                I = zeros(t+1, 1);
                for tau=1:t+1
                    Non(tau) = sum(cells_hist{tau});
                    I(tau) = moranI(cells_hist{tau}, a0*dist);
                end
                
                % save initial I
                %I_ini(k) = I(1);
                %p_out(k) = Non(end)/N;
                %eq_times1(k) = numel(Non) - 1;
                
                % Plot trajectory
                figure(h);
                hold on
                plot_sim = plot(Non/N, I, 'r-', 'LineWidth', 1.5);
                plot(Non(1)/N, I(1), 'ro', 'LineWidth', 2);
                plot(Non(end)/N, I(end), 'rx', 'LineWidth', 2);
                plot_sim.Color(4) = 1; % transparency
            end
        end
    end
end

%% New Langevin equation 
%{
% define gradient
delh_delp = @(p,I) (Con-1)*2*fN*(I-1).*(2*p - 1) -(Con+1)*(fN+1)+2*K;
delh_delI = @(p,I) 2*fN*(Con - 1)*p.*(p - 1);

% trajectory variables
maxsteps = 100;
p_stoch = zeros(nruns, maxsteps+1);
I_stoch = zeros(nruns, maxsteps+1);
pnoise = zeros(nruns, maxsteps+1);
Inoise = zeros(nruns, maxsteps+1);
eq_times2 = zeros(nruns, 1);

% If plotting multiple trajectories from same p_ini
%test = reshape(repmat(round(N*p_ini)/N, [10 1]), [30 1]);
%I_ini = 0;

pnoise(:,1) = 0;
Inoise(:,1) = 0; %normrnd(0, 0.02, [nruns, 1]); % take larger noise for initial spread
p_stoch(:,1) = p_ini' + pnoise(:,1);
%p_stoch(:,1) = test;
I_stoch(:,1) = I_ini + Inoise(:,1);

Peq_all = zeros(nruns, maxsteps);
for i=1:nruns
    for j=1:maxsteps 
        pt = p_stoch(i, j);
        It = I_stoch(i, j);
        [p_new, I_new, Peq, p_lim] = EOM_update_Langevin(pt, It, N, Con, K, fN, gN, delta, sigma, 0);
        Peq_all(i,j) = Peq;
        
        % system out of range?
        if p_lim == 1 || p_lim == -1
            p_stoch(i, j+1:end) = (p_lim+1)/2;
            I_stoch(i, j+1:end) = 0;
        end
        
        % system in equilibrium?
        if rand < Peq
            eq_times2(i) = j-1;
            p_stoch(i, j+1:end) = pt;
            I_stoch(i, j+1:end) = It;
            break
        else
            p_stoch(i, j+1) = p_new;
            I_stoch(i, j+1) = I_new;
        end
    end
end

% Plot Langevin trajectories
figure(h);
hold on
%h1a = plot(p_stoch', I_stoch', 'g-', 'LineWidth', 1.5);
handles = cell(numel(p_ini), 1);
for pl=1:nruns
    handles{pl} = plot(p_stoch(pl,:)', I_stoch(pl,:)', 'g-', 'LineWidth', 2);
    handles{pl}.Color(4) = 0.6;
end   

% plot start and end points
plot(p_stoch(:,1), I_stoch(:,1), 'go', 'LineWidth', 2);
plot(p_stoch(:,end), I_stoch(:,end), 'gx', 'LineWidth', 2);

%% set figure properties
set(gca, 'FontSize', 40);
xticks([0 1]);
yticks([0 1]);
set(c, 'XTick', [0 1]);
%set(gcf, 'CLim', [0 1])
%}
%% Save figure
qsave = 1;
if qsave
    save_path = 'H:\My Documents\Multicellular automaton\figures\random_positions';
    fname_str = strrep(sprintf('Peq_trajectories_N%d_a0_%.1f_K_%d_Con_%d_noise%.1f_%s', ...
        N, a0, K, Con, noise, labels{choice}), '.', 'p');
    fname = fullfile(save_path, fname_str);
    save_figure(h, 10, 8, fname, '.pdf');
end

%}