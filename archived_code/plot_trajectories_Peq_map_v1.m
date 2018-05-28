%% Plots a map of Peq in p, I space
clear all
close all
clc
%% Parameters
gridsize = 15;
N = gridsize^2;
a0 = 0.5;
Rcell = 0.2*a0;
K = 18; %12 20
Con = 11;

% load fN, gN
[dist, pos] = init_dist_hex(gridsize, gridsize);
dist_vec = dist(1,:);
r = a0*dist_vec(dist_vec>0); % exclude self influence
fN = sum(sinh(Rcell)*sum(exp(Rcell-r)./r)); % calculate signaling strengthN
gN = sum(sum((sinh(Rcell)*exp(Rcell-r)./r).^2)); % calculate signaling strength

% initial conditions
I0 = 0;
n1 = floor(0.1*N);
n_all = (n1:3:N-n1);
I_all = 0:0.05:0.3;
dI = 0.01;

%% Calculate P_eq
p = (0:N)/N;
I = linspace(-0.15, 1, 2*116);
[pm, Im] = meshgrid(p,I);

muon = fN*(Con.*pm + 1 - pm + (Con-1).*(1-pm).*Im);
muoff = fN*(Con.*pm + 1 - pm - (Con-1).*pm.*Im);

muon2 = fN*(Con.*pm + 1 - pm);
muoff2 = muon2;
%muon2 = Con*fN + (Con-1)/2*(2.*(1-pm).*Im + (pm-1)./pm*fN);
%muon2 = (Con+1)/2*fN + (Con-1)/2*((4*pm.*(1-pm).*Im + (2*pm-1)*fN)+(2*pm-1))./(2*pm);% 
%muoff2 = (Con+1)/2*fN - (Con-1).*pm.*Im;

kappa = sqrt((Con-1)^2*(gN)*pm.*(1-pm));
sigmaon = kappa;
sigmaoff = kappa;

zon = (K - Con - muon)./sigmaon;
zoff = (K - 1 - muoff)./sigmaoff;
zon2 = (K - Con - muon2)./sigmaon;
zoff2 = (K - 1 - muoff2)./sigmaoff;

Poffoff = normcdf(zoff);
Ponon = 1-normcdf(zon);
Poffoff2 = normcdf(zoff2);
Ponon2 = 1-normcdf(zon2);

pe = (Ponon.^pm .* Poffoff.^(1-pm)).^N;
pe2 = (Ponon2.^pm .* Poffoff2.^(1-pm)).^N;

%% Plot Ponon / Poffoff
%{
h=figure(2);
imagesc(p, I, Ponon);
title('Ponon');
set(gca, 'YDir', 'Normal');
colorbar;

h=figure(3);
imagesc(p, I, Poffoff);
title('Poffoff');
set(gca, 'YDir', 'Normal');
colorbar;
%}
%% Plot result
h=figure(1);
hold on
imagesc(p, I, pe);
set(gca, 'YDir', 'Normal');
c = colorbar;
xlabel('p');
ylabel('I');
ylabel(c, 'P_{eq}');
caxis([0 1]);
set(gca, 'FontSize', 24);
set(gcf, 'Units','Inches', 'Position', [0 0 10 8]);
xlim([0 1]);
ylim([I(1) I(end)]);

%%
h=figure(2);
hold on
imagesc(p, I, pe2);
set(gca, 'YDir', 'Normal');
c = colorbar;
xlabel('p');
ylabel('I');
ylabel(c, 'P_{eq}');
caxis([0 1]);
set(gca, 'FontSize', 24);
set(gcf, 'Units','Inches', 'Position', [0 0 10 8]);
xlim([0 1]);
ylim([I(1) I(end)]);
%% Plot gradient vector field
% define gradient
delh_delp = @(p,I) (Con-1)*2*fN*(I-1).*(2*p - 1) -(Con+1)*(fN+1)+2*K;
delh_delI = @(p,I) 2*fN*(Con - 1)*p.*(p - 1);

figure(h);
pv = (0:gridsize:N)/N;
Iv = -0.15:0.08:1;
[pm2, Im2] = meshgrid(pv, Iv);
vf_x = -delh_delp(pm2, Im2);
vf_y = -delh_delI(pm2, Im2);
quiver(pm2, Im2, vf_x, vf_y, 'LineWidth', 1, 'AutoScaleFactor', 1.2,...
    'Color', [0.5 0.5 0.5]);
%}
%% Load simulated trajectories

% Directory & filename of saved trajectories
foldername = ''; %strrep(sprintf('a0_%.1f_K%d_Con%d', a0, K, Con), '.', 'p');
path = fullfile('H:\My Documents\Multicellular automaton\data\dynamics\no_noise', foldername);
listing = dir(path);
num_files = numel(listing)-2; %first two entries are not useful

straux = '(\d+)';
%fpattern = sprintf('N%d_n%s_I0_%s_a0_%d_K_%d_Con_%d_neq_%s_t_%s-v%s', N, straux,...
%    straux, 10*a0, K, Con, straux, straux, straux);
fpattern = strrep(sprintf('N%d_n%s_I%.1f_neq_%s_a0_%.2f_K_%d_Con_%d_t_%s-v%s', N, straux,...
    I0, straux, a0, K, Con, straux, straux), '.', 'p');

% variables 
I_ini = [];
p_out = []; % store p_out values
eq_times1 = [];

nruns = 0;
names = {};
for i = 1:num_files
    filename = listing(i+2).name;
    % remove extension and do not include txt files
    [~,name,ext] = fileparts(filename);
    if strcmp(ext, '.mat')
        names{end+1} = name;
        [tokens, ~] = regexp(name, fpattern, 'tokens', 'match');
        if numel(tokens) > 0
            disp(name);
            nruns = nruns + 1;
            
            load(fullfile(path,name), 'Non', 'I', 'mom');
            iniON = Non(1);
            % save initial I
            %I_ini(k) = I(1);
            %p_out(k) = Non(end)/N;
            %eq_times1(k) = numel(Non) - 1;
            % Plot trajectory
            figure(h); %figure(1);
            hold on
            plot_sim = plot(Non/N, I, 'r-', 'LineWidth', 0.3);
            plot(Non(1)/N, I(1), 'ro', 'LineWidth', 1);
            plot(Non(end)/N, I(end), 'rx', 'LineWidth', 1);
            plot_sim.Color(4) = 1; % transparency
            %if count >= nruns % limit trajectory number
            %    break
            %end
        end
    end
end

%% Plot box with initial range
%figure(h);
%rectangle('Position', [n1/N I_all(1) (N-2*n1)/N I_all(end)-I_all(1)], 'EdgeColor', 'r')

%% Save figure
qsave = 0;
if qsave
    fname_str = sprintf('N%d_a0%d_K_%d_Con_%d_sweep_pI',...
        N, 10*a0, K, Con);
    i=0;
    fname = fullfile(pwd, 'figures', 'generate_I', 'trajectories',...
        strcat(fname_str,'_Peq_map_w_traj','-v',int2str(i)));
    while exist(fname, 'file') == 2
            i=i+1;
            fname = fullfile(pwd, 'figures', 'generate_I', 'trajectories',...
                strcat(fname_str,'_Peq_map_w_traj','-v',int2str(i)));
    end
    save_figure_pdf(h, 10, 8, fname);
    save_figure_eps(h, 10, 8, fname);
end
%}