%% Analyze trajectory output times and periods
clear all
close all
set(0, 'defaulttextinterpreter', 'latex');
%%
load_path = 'L:\BN\HY\Shared\Yiteng\two_signals\LHS_sample 1';
save_path_fig = 'H:\My Documents\Multicellular automaton\figures\two_signals\analyze_trajectories_LHS_sampling';

% Parameters
nruns = 100;
gz = 15;
N = gz^2;
M_int = [0 1; -1 0];
InitiateI = 0;
dist = init_dist_hex(gz, gz);
a0 = 1.5;
lambda = [1 1.2];
rcell = 0.2;
Rcell = rcell*a0;
dist_vec = a0*dist(1,:);
r = dist_vec(dist_vec>0); % exclude self influence
fN = zeros(2,1);
fN(1) = sum(sinh(Rcell)*sum(exp((Rcell-r)./lambda(1)).*(lambda(1)./r)) ); % calculate signaling strength
fN(2) = sum(sinh(Rcell)*sum(exp((Rcell-r)./lambda(2)).*(lambda(2)./r)) ); % calculate signaling strength

% Load LHS parameters
LHS_path = 'H:\My Documents\Multicellular automaton\data\two_signals\LHS_sample';
nvals = 30;
v = 1;
fname_str = sprintf('LHS_values_K12_K21_Con1_Con2_nvals%d-v%d',...
        nvals, v);
load(fullfile(LHS_path, fname_str), 'K12_all', 'K21_all', 'Con1_all', 'Con2_all');

% vars to store
t_out = zeros(nvals, nruns);
periods = zeros(nvals, nruns);
phase = zeros(1, nvals); % 0: oscillatory (possibly), 1: all ON/OFF phase
uniform = zeros(nvals, nruns); % final lattice uniform? 0: no, 1: yes

% loop over samples
for idx_s=1:nvals
    K = [0 K12_all(idx_s); K21_all(idx_s) 0];
    Con = [Con1_all(idx_s) Con2_all(idx_s)];
    
    % determine phase
    R12_ON = (1+fN(2) - K(1,2) ) > 0; % Everything ON
    R12_OFF = ((1+fN(2))*Con(2) - K(1,2) ) < 0; % Everything OFF
    R21_ON = (1+fN(1) - K(2,1) ) > 0; % Everything ON
    R21_OFF = ((1+fN(1))*Con(1) - K(2,1) ) < 0; % Everything OFF
    phase(idx_s) = (R12_ON | R12_OFF | R21_ON | R21_OFF);
    
    folder = '';
    listing = dir(fullfile(load_path, folder));
    num_files = numel(listing)-2; %first two entries are not useful

    d = '(\d+)';
    d1 = '0p(\d)0';
    d2 = '(\d+|Inf)';
    
    sim_ID = 'two_signal_mult';
    I_ini_str = '';
    pattern = strrep(sprintf(...
        '%s_N%d_initiateI%d%s_K12_%.2f_K21_%.2f_Con1_%.2f_Con2_%.2f_t_out_%s_period_%s',...
        sim_ID, N, InitiateI, I_ini_str, K(1,2), K(2,1), Con(1), Con(2),...
        d, d2), '.', 'p');
    disp(pattern);
    count = 0;
    for i = 1:num_files
        filename = listing(i+2).name;
        % remove extension and do not include txt files
        [~,name,ext] = fileparts(filename);
        if strcmp(ext, '.mat')
            [~, tokens] = regexp(name, pattern, 'match', 'tokens');
            if ~isempty(tokens)
                disp(name);
                load(fullfile(load_path, folder, strcat(name, ext)), 'cells_hist');
                cells = cells_hist{end};
                uniform(idx_s, count+1) = all(all(cells==0)+all(cells==1));
                t_out(idx_s, count+1) = str2double(tokens{1}{1});
                periods(idx_s, count+1) = str2double(tokens{1}{2});
                count = count + 1;
                %names{count} = name;
            end
        end
    end
end

%% (1) display LHS samples (as table?)


%% (2) Distribution of periods
f_periodic = sum(periods~=Inf, 2)/nruns;
f_period4 = sum(periods==4, 2)/nruns;

%h2=figure(2);
%hold on
%plot(1:nvals, f_periodic);
%plot(1:nvals, f_period4);

idx = (phase==0); % select only phase=0
% histogram
h21 = figure(21);
hold on
C1 = categorical(periods(idx,:), unique(periods(idx,:)));
C2 = categorical(periods(~idx,:), unique(periods(~idx,:)));
histogram(C1);
histogram(C2);
xlabel('Period');
ylabel('Count');
title(sprintf('Total %d runs', numel(periods)));
legend({'Normal phase', 'All ON/OFF'}, 'Location', 'eastoutside');

%% (3) Uniform lattice oscillations vs. individual oscillations

idx = (phase==0); % select only phase=0
C = categorical(uniform(idx,:), [0 1], {'Non-uniform', 'Uniform'});

h3 = figure(3);
histogram(C);
title(sprintf('Total %d runs', sum(idx)*nruns));

h31 = figure(31);
xdata = 1:sum(idx);
ydata = sum(uniform(idx,:), 2)/nruns;
plot(xdata, ydata, 'o-');
xlabel('Parameter set');
ylabel('Fraction uniform');
%% Plot together with phase diagram
% figuring out color map is a pain though... white=uniform,
% black=non-uniform?
% Also, each data point is a pairs of points in the two plots. The
% connections are not shown
cdata = sum(uniform(idx,:), 2)/nruns;
plot_phase_diagram(gz, a0, rcell, M_int, zeros(2), zeros(2,1), lambda(2));
%%
hold on
%h32 = figure(32);
subplot(2,1,1);
map = [0,0,0
    1,1,1];
colormap(map);
scatter(K21_all(idx), Con1_all(idx), 50, cdata)
text(K21_all(idx)+0.5, Con1_all(idx)+0.5, sprintfc('%d', 1:sum(idx)));
subplot(2,1,2);
scatter(K12_all(idx), Con2_all(idx), 50, cdata)
text(K12_all(idx)+0.5, Con2_all(idx)+0.5, sprintfc('%d', 1:sum(idx)));
