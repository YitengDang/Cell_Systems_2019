%% Analyze saved trajectories across a range of pin
clear all
close all
set(0, 'defaulttextinterpreter', 'latex');
%%
%path = 'H:\My Documents\Multicellular automaton\data\two_signals\time_evolution\III_chaotic_scan_p_ini_batch2';
path = 'D:\Multicellularity\data\two_signals\time_evolution\vs_pini_batch3';
listing = dir(path);
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

%% Load data
N = 225;
nruns = 80;
tmax = 10000;

p1 = 0:0.1:1;
p2 = 0:0.1:1;
N1 = round(p1*N);
N2 = round(p2*N);

filecount = zeros(numel(p1), numel(p2));
t_out_all = zeros(numel(p1), numel(p2), nruns); % final times
period = zeros(numel(p1), numel(p2), nruns); % periodicity test
t_onset = zeros(numel(p1), numel(p2), nruns); 

pattern = '.'; % '.' = anything 'p_ini\dp\d{2}_0p30';
for i=1:numel(names)
    if isempty(regexp(names{i}, pattern, 'once')) % only load files matching a certain pattern
        continue
    else
        disp(names{i});
        load(fullfile(path, names{i}));
        p_ini = save_consts_struct.p_ini;
        N_ini = round(p_ini*N);
        %disp(p_ini);
        idx1 = find(N_ini(1) == N1, 1);
        idx2 = find(N_ini(2) == N2, 1);
        if ~isempty(idx1)&&~isempty(idx2)
            filecount(idx1, idx2) = filecount(idx1, idx2) + 1;
            idx3 = filecount(idx1, idx2);
            t_out_all(idx1, idx2, idx3) = t_out;
            [period(idx1, idx2, idx3), t_onset(idx1, idx2, idx3)] = ...
                periodicity_test_v2(cells_hist); % periodicity test
        end
    end
end

%% Save the loaded data
fname_str = sprintf('III_chaotic_scan_p_ini_batch2_data_nruns%d', nruns);
save_path = 'H:\My Documents\Multicellular automaton\data\two_signals\time_evolution';
%save_vars = {'p1', 'p2', 'filecount', 't_out_all', 'period'};
save(fullfile(save_path, strcat(fname_str, '.mat') ), 'p1', 'p2', 'filecount', 't_out_all', 'period', 't_onset') %, save_vars);

% folder for saving figures
save_path_fig = 'H:\My Documents\Multicellular automaton\figures\two_signals\analyze_trajectories_vs_pin';

%% Load data
%fname_str = 'III_chaotic_scan_p_ini_data';
%load_path = 'H:\My Documents\Multicellular automaton\data\two_signals\time_evolution';
%load(fullfile(load_path, strcat(fname_str, '.mat')));
%% Analyze t_out
% (1) <t_out>
%t_out_mean = mean(t_out_all, 3);

% (2) # trajectories reaching tmax
t_max_reached = (t_out_all == tmax);
t_max_count = sum(t_max_reached, 3);
h2 = figure(2);
imagesc(p1, p2, t_max_count/nruns);
set(gca, 'YDir', 'Normal');
c = colorbar;
caxis([0 1]);
%colormap('parula');
xlabel('$$p_1$$')
ylabel('$$p_2$$')
ylabel(c, 'fraction');
set(gca, 'FontSize', 20);

qsave = 1;
if qsave
    fname = fullfile(save_path_fig, strcat(fname_str, '_frac_reaching_tmax'));
    save_figure(h2, 10, 8, fname, '.pdf');
end

%%
% (3) <t_out> over trajectories that don't reach tmax
t_out_mean_2 = zeros(numel(p1), numel(p2));
for i=1:numel(t_out_all)
    if ~t_max_reached(i)
        [i1,i2,i3] = ind2sub(size(t_out_all), i); 
        t_out_mean_2(i1, i2) = t_out_mean_2(i1, i2) + t_out_all(i1,i2,i3);
    end
end
t_out_mean_ = t_out_mean_2./(nruns - t_max_count);

h3 = figure(3);
imagesc(p1, p2, t_out_mean_2);
set(gca, 'YDir', 'Normal');
c = colorbar;
xlabel('$$p_1$$')
ylabel('$$p_2$$')
ylabel(c, '$$\langle t_{out} |  t_{out} < t_{max} \rangle$$', ...
    'Interpreter', 'latex', 'FontSize', 20);
set(gca, 'FontSize', 20);
caxis([0 tmax]);
%ylim(c, [0 tmax]);

qsave = 0;
if qsave
    fname = fullfile(save_path_fig, strcat(fname_str, '_mean_t_out_subset'));
    save_figure(h3, 10, 8, fname, '.pdf');
end
%% Analyze periods
% number of chaotic trajectories
n_chaotic = sum((period(:)==Inf).*(t_out_all(:) == tmax));

% data
period_data = period(period~=0);
uniq = unique(period_data);
C = categorical(period_data, uniq);

n_periodic = sum(period_data(:)~=Inf);

% distribution of periods
h4 = figure(4);
histogram(C);
xlabel('Period');
ylabel('Count');
title(sprintf('%d trajectories, %d periodic, %d chaotic',...
    sum(sum(filecount)), n_periodic, n_chaotic));
set(gca, 'FontSize', 20);

qsave = 1;
if qsave
    fname = fullfile(save_path_fig, strcat(fname_str, '_periods_hist_chaotic'));
    save_figure(h4, 10, 8, fname, '.pdf');
end


%% average period (count only found periods)

period_mean = zeros(numel(p1), numel(p2));
period_count = zeros(numel(p1), numel(p2));
for i=1:numel(t_out_all)
    if period(i)>0 && period(i)<Inf
        [i1,i2,i3] = ind2sub(size(period), i); 
        period_count(i1, i2) = period_count(i1, i2) + 1;
        period_mean(i1, i2) = period_mean(i1, i2) + period(i1,i2,i3);
    end
end
period_mean = period_mean./period_count;

h5 = figure(5);
imagesc(p1, p2, period_mean);
set(gca, 'YDir', 'Normal');
c = colorbar;
caxis([0 max(uniq(uniq<Inf))])
xlabel('$$p_1$$')
ylabel('$$p_2$$')
ylabel(c, 'Mean period (time steps)', ...
    'Interpreter', 'latex', 'FontSize', 20);
set(gca, 'FontSize', 20);

qsave = 0;
if qsave
    fname = fullfile(save_path_fig, strcat(fname_str, '_mean_period'));
    save_figure(h5, 10, 8, fname, '.pdf');
end