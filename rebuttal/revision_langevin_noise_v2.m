% Simulates Langevin equation based on gradient of Hamiltonian with random
% Gaussian noise
% Overlay simulated trajectories from the cellular automaton
clear all
close all
clc

%% Parameters
gridsize = 6;
N = gridsize^2;
a0 = 1.5;
Rcell = 0.2*a0;
%Klist = [10 10 14 19 16 18]; %(a0=1.5) [3 6 10 15 17 20]; % (a0=0.5) [10 10 14 19 16 18]; %[3]; 
%Conlist = [5 21 16 14 8 6]; %(a0=1.5) [24 21 21 20 14 14]; % (a0=0.5) [5 21 16 14 8 6]; %[24];
Klist = [6 20 12]; 
Conlist = [15 15 15];

%for idx = 1:numel(Klist)
idx = 1;
h = figure(idx);
hold on
%K = Klist(idx);
%Con = Conlist(idx);
K = 3;
Con = 24;

nruns = 30; % number of simulations
maxsteps = 1000; % number of steps per simulation
delta = 1/N; % % step size in Langevin equation, suggested 0.001
noise = 1; % 2; % Lagevin noise amplitude, suggested 2
alpha = 0; % sensing noise

% initial conditions
%p_ini = 0.65;
I_ini = 0;
p_ini = round( linspace(0,N, nruns) )/N;

% load fN
[dist, pos] = init_dist_hex(gridsize, gridsize);
dist_vec = dist(1,:);
r = a0*dist_vec(dist_vec>0); % exclude self influence
fN = sum(sinh(Rcell)*sum(exp(Rcell-r)./r)); % calculate signaling strengthN
gN = sum(sum((sinh(Rcell)*exp(Rcell-r)./r).^2)); % calculate signaling strength
%% Experiment with modifying drift in I
Ns = 6^2; %optimal N
fmax = 10; %max amplification
a = Ns*log(fmax);
Nrange = 25:400;
f_all = fmax.*exp(-a./Nrange);
%figure();
%plot(Nrange, f_all);
f = fmax.*exp(-a./N);

%% Overlay simulated trajectories
% load all files from a folder
%path = fullfile(pwd, 'rebuttal', 'trajectories', 'N121_pini_0p65_a0_0p5_K_16_Con_8');
path = fullfile(pwd, 'rebuttal', 'trajectories', 'N121_a0_1p5_figS4', sprintf('K%d_Con%d', K, Con));
%path = fullfile(pwd, 'rebuttal', 'trajectories', 'Fig_3_PRL_manuscript', sprintf('K%d_Con%d', K, Con));
listing = dir(path);
num_files = numel(listing)-2; %first two entries are not useful
count = 0;
% variables 
I_ini = zeros(numel(p_ini), 1);
p_out = zeros(numel(p_ini), 1); % store p_out values
eq_times1 = zeros(numel(p_ini), 1);
straux = '(\d+)';
fpattern = sprintf('N%d_n%s_neq_%s_a0%d_K_%d_Son_%d_t_%s-v%s',...
    N, straux, straux, 10*a0, K, Con, straux, straux);
for i = 1:num_files
    filename = listing(i+2).name;
    % remove extension and do not include txt files
    [~,name,ext] = fileparts(filename);
    if strcmp(ext, '.mat')
        [tokens, ~] = regexp(name, fpattern, 'tokens', 'match');
        if numel(tokens) > 0
            this_p_ini = str2num(tokens{1}{1})/N;
            if find(this_p_ini == p_ini)
                disp(name);
                count = count + 1;
                names{count} = name;
                load(fullfile(path,name), 'Non', 'I', 'mom');
                iniON = Non(1);
                % save initial I
                k = find(iniON/N == p_ini);
                I_ini(k) = I(1);
                p_out(k) = Non(end)/N;
                eq_times1(k) = numel(Non) - 1;
                % Plot trajectory
                h; %figure(1);
                hold on
                plot_sim = plot(Non/N, I, 'r-', 'LineWidth', 1.5);
                plot(Non(1)/N, I(1), 'ro');
                plot(Non(end)/N, I(end), 'rx');
                plot_sim.Color(4) = 0.6; % transparency
                %if count >= nruns % limit trajectory number
                %    break
                %end
            end
        end
    end
end

%figure(2);
%xlim([0 10]);
%% Simulation
% define gradient
hamiltonian = @(p, I) -0.5*(Con-1)*(1 + 4*fN.*p.*(1-p).*I + fN*(2*p-1).^2) ...
    -(2*p-1).*(0.5*(Con+1)*(1+fN) - K);
delh_delp = @(p,I) (Con-1)*2*fN*(I-1).*(2*p - 1) -(Con+1)*(fN+1)+2*K;
delh_delI = @(p,I) 2*fN*(Con - 1)*p.*(p - 1);

% trajectory variables
p_stoch = zeros(nruns, maxsteps+1);
I_stoch = zeros(nruns, maxsteps+1);
pnoise = zeros(nruns, maxsteps+1);
Inoise = zeros(nruns, maxsteps+1);
eq_times2 = zeros(nruns, 1);

pnoise(:,1) = 0;
Inoise(:,1) = 0; %normrnd(0, 0.02, [nruns, 1]); % take larger noise for initial spread
p_stoch(:,1) = p_ini' + pnoise(:,1);
I_stoch(:,1) = I_ini + Inoise(:,1);

% run simulation
%{
for i=1:nruns
    for j=1:maxsteps
        pt = p_stoch(i,j); 
        It = I_stoch(i,j);
        [peq, p_lim] = revision_update_Peq(pt, It, N, Con, K, fN, gN, alpha); %
        
        % system out of range?
        if p_lim == 1 || p_lim == -1 
            eq_times2(i) = j-1;
            p_stoch(i, j:end) = (p_lim+1)/2;
            I_stoch(i, j:end) = 0;
            break
        end
        
        % system in equilibrium?
        if rand() < peq 
            eq_times2(i) = j-1;
            p_stoch(i, j+1:end) = pt;
            I_stoch(i, j+1:end) = It;
            break
        else
            pnoise(i,j+1) = normrnd(0,noise)*delta;
            Inoise(i,j+1) = normrnd(0,noise)*delta;
            p_stoch(i,j+1) = p_stoch(i,j) - delh_delp(pt, It)*delta + pnoise(i,j+1);
            I_stoch(i,j+1) = I_stoch(i,j) - delh_delI(pt, It)*delta + Inoise(i,j+1);
        end
    end
end
%}
%
Peq_all = zeros(nruns, maxsteps);
for i=1:nruns
    for j=1:maxsteps 
        pt = p_stoch(i, j);
        It = I_stoch(i, j);
        [p_new, I_new, Peq, p_lim] = revision_update_Langevin(pt, It, N, Con, K, fN, gN, delta, noise, alpha);
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

%}
% Save trajectories
qsave = 1;
if qsave
    %fname_str = strrep(sprintf('N%d_pini_range_nruns%d_a0_%.1f_K_%d_Con_%d_delta_%.3f_noise_%d', ...
    %    N, nruns, a0, K, Con, delta, noise), '.', 'p');
    fname_str = strrep(sprintf('N%d_pini_range_nruns%d_a0_%.1f_K_%d_Con_%d_delta_%.3f_noise_%d_Idrift_f%.2f', ...
    N, nruns, a0, K, Con, delta, noise, f), '.', 'p');
    fname = fullfile(pwd, 'rebuttal', 'h_vector_field', 'Langevin_v2', fname_str);
    
    save(fname);
end

%% Plot result
%h=figure(1);
%hold on

% Plot vector field
pv = (0:gridsize:N)/N;
Iv = -0.1:0.1:1;
[pm, Im] = meshgrid(pv, Iv);
vf_x = -delh_delp(pm, Im);
vf_y = -delh_delI(pm, Im);
quiver(pm, Im, vf_x, vf_y, 'LineWidth', 1, 'AutoScaleFactor', 1.2);

% plot Langevin trajectories
h;
%h1a = plot(p_stoch', I_stoch', 'g-', 'LineWidth', 1.5);
handles = cell(numel(p_ini), 1);
for pl=1:numel(p_ini)
    handles{pl} = plot(p_stoch(pl,:)', I_stoch(pl,:)', 'g-', 'LineWidth', 1.5);
    handles{pl}.Color(4) = 0.5;
end   

% plot start and end points
plot(p_stoch(:,1), I_stoch(:,1), 'go', 'LineWidth', 1.5);
plot(p_stoch(:,end), I_stoch(:,end), 'gx', 'LineWidth', 1.5);

% plot layout
xlim([-0.2 1.2]);
ylim([-0.15 1.05]);
set(gca, 'FontSize', 24);
xlabel('p');
ylabel('I');

%% Save figure
qsave = 1;
if qsave
    %fname_str = strrep(sprintf('N%d_pini_%.2f_a0_%.1f_K_%d_Con_%d', ...
    %    N, p_ini, a0, K, Con), '.', 'p');
    %fname_str = strrep(sprintf('N%d_pini_range_nruns%d_a0_%.1f_K_%d_Con_%d_delta_%.3f_noise_%d', ...
    %    N, nruns, a0, K, Con, delta, noise), '.', 'p');
    fname_str = strrep(sprintf('N%d_pini_range_nruns%d_a0_%.1f_K_%d_Con_%d_delta_%.3f_noise_%d_Idrift_f%.2f', ...
    N, nruns, a0, K, Con, delta, noise, f), '.', 'p');
    fname = fullfile(pwd, 'rebuttal', 'h_vector_field', 'Langevin_v2', fname_str);
    save_figure_pdf(h, 10, 8, fname);
end

%end %----end K loop------

%% Plot locations of final values of p
figure();
hold on
plot(p_ini, p_out, 'r-');
plot(p_ini, p_stoch(:,end), 'g-');
xlabel('p_{in}');
ylabel( 'p_{out}');
%% Plot equilibration times
figure();
hold on 
plot(p_ini, eq_times1, 'r-');
plot(p_ini, eq_times2, 'g-');
xlabel('p_{in}');
ylabel('t_{eq}');
%% Plot h against t for Langevin trajectories
%{
h2=figure(2);
hold on
traj = 1:nruns; % trajectories to plot
plot(0:nsteps, hamiltonian(p_stoch(traj, :), I_stoch(traj, :) ) );
xlabel('t');
ylabel('h');
%}