% Simulates Langevin equation based on gradient of Hamiltonian with random
% Gaussian noise
% Overlay simulated trajectories from the cellular automaton
% v3: estimate noise parameters from Eduardo's equation, no free noise
% parameters
clear all
close all
clc
%% Parameters
gridsize = 15;
N = gridsize^2;
a0 = 0.5;
Rcell = 0.2*a0;
%Klist = [10 10 14 19 16 18]; %(a0=1.5) [3 6 10 15 17 20]; % (a0=0.5) [10 10 14 19 16 18]; %[3]; 
%Conlist = [5 21 16 14 8 6]; %(a0=1.5) [24 21 21 20 14 14]; % (a0=0.5) [5 21 16 14 8 6]; %[24];

%for idx = 1:numel(Klist)
%K = Klist(idx);
%Con = Conlist(idx);
K = 15;
Con = 8;
idx = 1;
h = figure(idx);
hold on

nruns = 30; % number of simulations
maxsteps = 1000; % number of steps per simulation
%delta = 0.01; % % step size in Langevin equation, suggested 0.001
%noise = 2; % 2; % Langevin noise amplitude, suggested 2
alpha = 0; % sensing noise

% initial conditions
p_ini = 0.5;
I_ini = 0;
%p_ini = round( linspace(0,N, nruns) )/N;

% load fN
[dist, pos] = init_dist_hex(gridsize, gridsize);
dist_vec = dist(1,:);
r = a0*dist_vec(dist_vec>0); % exclude self influence
fN = sum(Rcell*sum(exp(Rcell-r)./r)); % calculate signaling strengthN
gN = sum(sum((sinh(Rcell)*exp(Rcell-r)./r).^2)); % calculate signaling strength

%% Plot simulated trajectories
% load all files from a folder
%path = fullfile(pwd, 'rebuttal', 'trajectories', 'N121_pini_0p8_a0_0p5_K_16_Con_8');
path = fullfile(pwd, 'rebuttal\h_vector_field\data\a0_0p5_figS5', sprintf('K%d_Con%d', K, Con));

listing = dir(path);
num_files = numel(listing)-2; %first two entries are not useful
count = 0;
%I_ini = zeros(numel(p_ini), 1); % p_ini range
I_ini = zeros(nruns, 1); % p_ini fixed
for i = 1:num_files
    filename = listing(i+2).name;
    % remove extension and do not include txt files
    [~,name,ext] = fileparts(filename);
    if strcmp(ext, '.mat')
        count = count + 1;
        names{count} = name;
        load(fullfile(path,name), 'Non', 'I', 'mom');
        iniON = Non(1);
        % (1) fixed p_ini
        I_ini(count) = I(1);
        % comment out if loop below
        
        % (2) range of p_ini
        %del_n = round((1/nruns*N)); % take steps of del_n
        %if mod(iniON, del_n) == 0 % select a subset of trajectories
        %if find(iniON == round( linspace(0,N, nruns) ))
            % save initial I
            %k = find(iniON == round( linspace(0,N, nruns) ));
            %I_ini(k) = I(1);
            
            % Plot trajectory
            h;
            hold on
            plot_sim = plot(Non/N, I, 'r-', 'LineWidth', 1.5);
            plot(Non(1)/N, I(1), 'ro');
            plot(Non(end)/N, I(end), 'rx');
            plot_sim.Color(4) = 0.5; % transparency
            if count >= nruns % limit trajectory number
                break
            end
        %end
        % Plot hamiltonian trajectories
        %figure(2);
        %plot(0:numel(mom)-1, mom/N, 'b-');
    end
end

%figure(2);
%xlim([0 10]);

%% Simulation
% define gradient
%hamiltonian = @(p, I) -0.5*(Con-1)*(1 + 4*fN.*p.*(1-p).*I + fN*(2*p-1).^2) ...
%    -(2*p-1).*(0.5*(Con+1)*(1+fN) - K);
delh_delp = @(p,I) (Con-1)*2*fN*(I-1).*(2*p - 1) -(Con+1)*(fN+1)+2*K;
delh_delI = @(p,I) 2*fN*(Con - 1)*p.*(p - 1);

% trajectory variables
p_stoch = zeros(nruns, maxsteps+1);
I_stoch = zeros(nruns, maxsteps+1);
pnoise = zeros(nruns, maxsteps+1);
Inoise = zeros(nruns, maxsteps+1);

pnoise(:,1) = 0;
Inoise(:,1) = 0; %normrnd(0, 0.02, [nruns, 1]); % take larger noise for initial spread
p_stoch(:,1) = round(p_ini'*N)/N + pnoise(:,1);
I_stoch(:,1) = I_ini + Inoise(:,1);

% run simulation
%{
eq_times2 = zeros(nruns, 1);
for i=1:nruns
    for j=1:maxsteps
        pt = p_stoch(i,j); 
        It = I_stoch(i,j);

        [peq, p_lim] = revision_update_Peq(pt, It, N, Con, K, fN, gN, alpha); %
        if rand() < peq % system in equilibrium?
            if p_lim == 1 || p_lim == -1 %if system goes out of range
                p_stoch(i, j:end) = (p_lim+1)/2;
                I_stoch(i, j:end) = 0;
                eq_times2(i) = j-1; 
            else
                p_stoch(i, j:end) = pt;
                I_stoch(i, j:end) = It;
                eq_times2(i) = j-1; 
            end
            break
        end
        %--If not in equilibrium, update--
        % Calculate <Delta p> according to Eduardo's formula 
        muon = fN*(Con*pt + 1 - pt + (Con-1)*(1-pt)*It);
        muoff = fN*(Con*pt + 1 - pt - (Con-1)*pt*It);
        kappa = sqrt((Con-1)^2*(gN)*pt*(1-pt));
        alpha = 0; % noise in K
        sigmaon = sqrt(kappa^2+alpha^2);
        sigmaoff = sqrt(kappa^2+alpha^2);
        zon = (K - Con - muon)/sigmaon;
        zoff = (K - 1 - muoff)/sigmaoff;
        Poffoff = normcdf(zoff);
        Ponon = 1-normcdf(zon);
        delp = (1-pt)*(1-Poffoff) - pt*(1-Ponon);
        % match delta so that change in p is equal to <Delta p>
        delta = -delp/delh_delp(pt, It);
        
        % match noise with std(Delta p)
        noise = sqrt(pt*(1-Poffoff)*Poffoff + pt*(1-Ponon)*Ponon)/sqrt(N);
        %noise = 0.02; 
        
        % Update p, I
        pnoise(i,j+1) = normrnd(0,noise);
        Inoise(i,j+1) = normrnd(0,noise)*delh_delI(pt, It)/delh_delp(pt, It);
        p_stoch(i,j+1) = p_stoch(i,j) - delh_delp(pt, It)*delta + pnoise(i,j+1);
        I_stoch(i,j+1) = I_stoch(i,j) - delh_delI(pt, It)*delta + Inoise(i,j+1);
    end
end
%}
%
eq_times2 = zeros(nruns, 1);
Peq_all = zeros(nruns, maxsteps);
for i=1:nruns
    for j=1:maxsteps 
        pt = p_stoch(i, j);
        It = I_stoch(i, j);
        [p_new, I_new, Peq, p_lim] = revision_update_Langevin_no_noise_params(pt, It, N, Con, K, fN, gN, alpha);
        Peq_all(i,j) = Peq;
        
        % system out of range?
        if p_lim == 1 || p_lim == -1
            eq_times2(i) = j-1;
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
%% Plot Langevin result in same graph
%h=figure(1);
%hold on

% Plot vector field
h;
hold on
pv = (0:5:N)/N;
Iv = -0.1:0.1:1;
[pm, Im] = meshgrid(pv, Iv);
vf_x = -delh_delp(pm, Im);
vf_y = -delh_delI(pm, Im);
quiver(pm, Im, vf_x, vf_y, 'LineWidth', 1, 'AutoScaleFactor', 1.2);

%---plot Langevin trajectories---
% (1) All same p_ini
%h1a = plot(p_stoch', I_stoch', 'g-', 'LineWidth', 1.5);

% (2) range of different p_ini
handles = cell(nruns, 1);
for pl=1:nruns
    handles{pl} = plot(p_stoch(pl,:)', I_stoch(pl,:)', 'g-', 'LineWidth', 1.5);
    handles{pl}.Color(4) = 0.5;
end      
% --------------------------

% plot start and end points
h2a = plot(p_stoch(:,1), I_stoch(:,1), 'go', 'LineWidth', 1.5);
h2b = plot(p_stoch(:,end), I_stoch(:,end), 'gx', 'LineWidth', 1.5);
h2a.Color(4) = 0.5;
h2b.Color(4) = 0.5;

% plot layout
xlim([-0.2 1.2]);
ylim([-0.15 1.05]);
set(gca, 'FontSize', 24);
xlabel('p');
ylabel('I');

%% Save figure
qsave = 0;
if qsave
    fname_str = strrep(sprintf('N%d_pini_%.2f_a0_%.1f_K_%d_Con_%d', ...
        N, p_ini, a0, K, Con), '.', 'p');
    %fname_str = strrep(sprintf('N%d_pini_range_nruns%d_a0_%.1f_K_%d_Con_%d', ...
    %    N, nruns, a0, K, Con), '.', 'p');
    fname = fullfile(pwd, 'rebuttal', 'h_vector_field', 'Langevin_v3', fname_str);
    save_figure_pdf(h, 10, 8, fname);
end

%end %----end Klist loop---

%% Changes in p
%
%traj = 2;
figure();
hold on
for i=1:nruns
    plot(0:eq_times2(i), p_stoch(i, 1:eq_times2(i)+1) );
    plot(eq_times2(i), p_stoch(i, eq_times2(i)+1), 'x');
    %plot(0, p_ini(i), 'o');
end
%%
figure();
plot(p_ini, eq_times2, 'o-');
%}
%% Plot h against t for Langevin trajectories
%{
h2=figure(2);
hold on
traj = 1:nruns; % trajectories to plot
plot(0:nsteps, hamiltonian(p_stoch(traj, :), I_stoch(traj, :) ) );
xlabel('t');
ylabel('h');
%}