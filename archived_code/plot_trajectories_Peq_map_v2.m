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
p0 = 0.5;
iniON = round(p0*N);
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
%% Plot Peq map
h1=figure(1);
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

%% Plot gradient vector field
% define gradient
delh_delp = @(p,I) (Con-1)*2*fN*(I-1).*(2*p - 1) -(Con+1)*(fN+1)+2*K;
delh_delI = @(p,I) 2*fN*(Con - 1)*p.*(p - 1);

figure(h1);
pv = (0:gridsize:N)/N;
Iv = -0.15:0.08:1;
[pm2, Im2] = meshgrid(pv, Iv);
vf_x = -delh_delp(pm2, Im2);
vf_y = -delh_delI(pm2, Im2);
quiver(pm2, Im2, vf_x, vf_y, 'LineWidth', 1, 'AutoScaleFactor', 1.2,...
    'Color', [0.5 0.5 0.5]);
%}
%% Load simulated trajectories

% Path to search for the saved data. It searchs by the name, defined by the
% parameters chosen
path_CA = 'H:\My Documents\Multicellular automaton\data\dynamics\no_noise\batch1';
straux = '(\d+)';

% filename pattern (regexp notation)
%fpattern_CA = sprintf('N%d_n%d_neq_%s_a0%d_K_%d_Son_%d_t_%s-v%s',...
%    N, iniON, straux, 10*a0, K, Con, straux, straux);
%fpattern_CA = sprintf('N%d_n%s_neq_%s_a0%d_K_%d_Son_%d_t_%s-v%s',...
%    N, straux, straux, 10*a0, K, Con, straux, straux);
%fpattern_CA = strrep(sprintf('N%d_n%s_I%.1f_neq_%s_a0_%.2f_K_%d_Con_%d_t_%s', ...
%    N, straux, I0, straux, a0, K, Con, straux), '.', 'p');

% Get all file names in the directory
listing = dir(path_CA);
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

% Initialize all variables and figures
pini = zeros(numel(names),1);
pend = pini;
Iini = pini;
Iend = pini;
hini = pini;
hend = pini;
tend = pini;
for i = 1:numel(names)
    % first get the filename given by the student
    [tokens, ~] = regexp(names{i},fpattern_CA,'tokens','match');
    if numel(tokens) > 0
        disp(names{i}) % displays the file name
        % load the data
        load(fullfile(path_CA, strcat(names{i},'.mat')), 'Non', 'fN', 'I', 't');
        % Get the sequence of p
        p = Non/N;
        
        figure(h1)
        % Plot the line
        plot(p, I, 'Color', 'r')

        % Save the initial and end points for plotting the marks later
        pend(i) = p(end);
        Iend(i) = I(end);
        pini(i) = p(1);
        Iini(i) = I(1);
        tend(i) = t;
    end
end

%%
% Path to search for the saved data. It searchs by the name, defined by the
% parameters chosen
path_MC = 'H:\My Documents\Multicellular automaton\data\dynamics\no_noise_MC\dp_stochastic\fixed_pin';
straux = '(\d+)';
%fpattern_MC = sprintf('N%d_n%d_neq_%s_a0%d_K_%d_Son_%d_t_%s_montecarlo_v3-v%s', ...
%    N, iniON, straux, 10*a0, K, Con, straux, straux);
fpattern_MC = strrep(sprintf('N%d_n%d_I0_%.1f_neq_%s_a0_%.2f_K_%d_Con_%d_t_%s_montecarlo_v3_dp_stoch', ...
    N, iniON, I0, straux, a0, K, Con, straux), '.', 'p');

% Get all file names in the directory
listing = dir(path_MC);
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

% Initialize all variables and figures
pini = zeros(numel(names),1);
pend = pini;
Iini = pini;
Iend = pini;
Iini_2 = pini;
Iend_2 = pini;
hini = pini;
hend = pini;
tend = pini;
for i = 1:numel(names)
    % first get the filename given by the student
    [tokens, ~] = regexp(names{i},fpattern_MC,'tokens','match');
    if numel(tokens) > 0
        disp(names{i}) % displays the file name
        
        % load the data
        load(fullfile(path_MC,strcat(names{i},'.mat')), 'Non', 'fN', 'I', 'I_2', 't');
        % Get the sequence of p
        p = Non/N;
        figure(h1)
        % Plot the line
        %plot(p, I, 'Color', 'r')
        plot(p, I_2, 'Color', 'b'); % blue: MC simulations (v3)

        % Save the initial and end points for plotting the marks later
        pend(i) = p(end);
        Iend(i) = I(end);
        Iend_2(i) = I_2(end);
        pini(i) = p(1);
        Iini(i) = I(1);
        Iini_2(i) = I_2(1);
        tend(i) = t;
    end
end

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