%% Plots a map of Peq in p, I space
% Plot trajectories (simulated, Langevin) on top of it
clear all
close all
clc

%% Parameters
gridsize = 11;
N = gridsize^2;
a0 = 0.5;
Rcell = 0.2*a0;
K = 16;
Con = 8;
%Klist =  [3 6 10 15 17 20]; %(a0=1.5) [3 6 10 15 17 20]; % (a0=0.5) [10 10 14 19 16 18]; %[3]; 
%Conlist = [24 21 21 20 14 14]; %(a0=1.5) [24 21 21 20 14 14]; % (a0=0.5) [5 21 16 14 8 6]; %[24];

% load fN, gN
[dist, pos] = init_dist_hex(gridsize, gridsize);
dist_vec = dist(1,:);
r = a0*dist_vec(dist_vec>0); % exclude self influence
fN = sum(Rcell*sum(exp(Rcell-r)./r)); % calculate signaling strengthN
gN = sum(sum((sinh(Rcell)*exp(Rcell-r)./r).^2)); % calculate signaling strength

% --------For loop K, Con--------------
%for i=1:numel(Klist)
%    K = Klist(i);
%    Con = Conlist(i);
%% Calculate P_eq
p = (0:N)/N;
I = linspace(-0.15, 1, 2*116);
[pm, Im] = meshgrid(p,I);

muon = fN*(Con.*pm + 1 - pm + (Con-1).*(1-pm).*Im);
muoff = fN*(Con.*pm + 1 - pm - (Con-1).*pm.*Im);
kappa = sqrt((Con-1)^2*(gN)*pm.*(1-pm));
sigmaon = kappa;
sigmaoff = kappa;

zon = (K - Con - muon)./sigmaon;
zoff = (K - 1 - muoff)./sigmaoff;

Poffoff = normcdf(zoff);
Ponon = 1-normcdf(zon);

pe = (Ponon.^pm .* Poffoff.^(1-pm)).^N;
%% Plot result
h=figure(1);
imagesc(p, I, pe);
set(gca, 'YDir', 'Normal');
c = colorbar;
xlabel('p');
ylabel('I');
ylabel(c, 'P_{eq}');
caxis([0 1]);
set(gca, 'FontSize', 24);
set(gcf, 'Units','Inches', 'Position', [0 0 10 8]);
%% Load simulated trajectories
% p values to take
%p_ini = (0:5:N)/N;
p_ini = [0.5 0.65 0.8];
nruns = 0;

% load all files from a folder
%path = fullfile(pwd, 'rebuttal', 'trajectories', 'N121_pini_0p65_a0_0p5_K_16_Con_8');
%path = fullfile(pwd, 'rebuttal\h_vector_field\data\a0_1p5_figS4', sprintf('K%d_Con%d', K, Con));
%path = fullfile(pwd, 'rebuttal', 'trajectories', 'Fig_3_PRL_manuscript', sprintf('K%d_Con%d', K, Con));
path = fullfile(pwd, 'rebuttal', 'trajectories', 'Peq_figure');
listing = dir(path);
num_files = numel(listing)-2; %first two entries are not useful
count = 0;
%I_ini = zeros(numel(p_ini), 1);
%p_out = zeros(numel(p_ini), 1); % store p_out values
eq_times1 = zeros(numel(p_ini), 1);
names = {};
for i = 1:num_files
    filename = listing(i+2).name;
    % remove extension and do not include txt files
    [~,name,ext] = fileparts(filename);
    if strcmp(ext, '.mat')
        count = count + 1;
        names{count} = name;
        load(fullfile(path,name), 'Non', 'I', 'mom');
        iniON = Non(1);
        if find(iniON == round(N*p_ini))
            % save initial I
            k = find(iniON == round(N*p_ini) );
            %I_ini(k) = I(1);
            %p_out(k) = Non(end)/N;
            eq_times1(k) = numel(Non) - 1;
            nruns = nruns + 1;
            % Plot trajectory
            h; %figure(1);
            hold on
            plot_sim = plot(Non/N, I, 'r-', 'LineWidth', 2);
            plot(Non(1)/N, I(1), 'ro');
            plot(Non(end)/N, I(end), 'rx');
            plot_sim.Color(4) = 1; % transparency
            %if count >= nruns % limit trajectory number
            %    break
            %end
        end
    end
end

%% Plot gradient vector field
% define gradient
delh_delp = @(p,I) (Con-1)*2*fN*(I-1).*(2*p - 1) -(Con+1)*(fN+1)+2*K;
delh_delI = @(p,I) 2*fN*(Con - 1)*p.*(p - 1);

h;
pv = (0:gridsize:N)/N;
Iv = -0.15:0.08:1;
[pm2, Im2] = meshgrid(pv, Iv);
vf_x = -delh_delp(pm2, Im2);
vf_y = -delh_delI(pm2, Im2);
quiver(pm2, Im2, vf_x, vf_y, 'LineWidth', 1, 'AutoScaleFactor', 1.2,...
    'Color', [0.5 0.5 0.5]);

%% New Langevin equation 
%
test = reshape(repmat(round(N*p_ini)/N, [10 1]), [30 1]);
I_ini = 0;

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

pnoise(:,1) = 0;
Inoise(:,1) = normrnd(0, 0.02, [nruns, 1]); %normrnd(0, 0.02, [nruns, 1]); % take larger noise for initial spread
%p_stoch(:,1) = p_ini' + pnoise(:,1);
p_stoch(:,1) = test;
I_stoch(:,1) = I_ini + Inoise(:,1);

Peq_all = zeros(nruns, maxsteps);
for i=1:nruns
    for j=1:maxsteps 
        pt = p_stoch(i, j);
        It = I_stoch(i, j);
        [p_new, I_new, Peq, p_lim] = revision_update_Langevin_no_noise_params(pt, It, N, Con, K, fN, gN, 0);
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
h;
%h1a = plot(p_stoch', I_stoch', 'g-', 'LineWidth', 1.5);
handles = cell(numel(p_ini), 1);
for pl=1:nruns
    handles{pl} = plot(p_stoch(pl,:)', I_stoch(pl,:)', 'g-', 'LineWidth', 1.5);
    handles{pl}.Color(4) = 0.6;
end   

% plot start and end points
plot(p_stoch(:,1), I_stoch(:,1), 'go', 'LineWidth', 1.5);
plot(p_stoch(:,end), I_stoch(:,end), 'gx', 'LineWidth', 1.5);
%}
%% Save figure
qsave = 0;
if qsave
    %fname_str = strrep(sprintf('N%d_pini_%.2f_a0_%.1f_K_%d_Con_%d', ...
    %    N, p_ini, a0, K, Con), '.', 'p');
    fname_str = strrep(sprintf('N%d_pini_range_nruns%d_a0_%.1f_K_%d_Con_%d_zero_param', ...
        N, nruns, a0, K, Con), '.', 'p');
    fname = fullfile(pwd, 'rebuttal', 'Peq_map_w_trajectories', strcat('Peq_vf_traj_', fname_str));
    save_figure_pdf(h, 10, 8, fname);
    save_figure_eps(h, 10, 8, fname);
end

%-------end for loop K, Con---------
%end