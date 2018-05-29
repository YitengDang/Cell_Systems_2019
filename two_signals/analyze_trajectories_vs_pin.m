%% Analyze saved trajectories across a range of pin
clear all
close all
set(0, 'defaulttextinterpreter', 'latex');
%%
path = 'H:\My Documents\Multicellular automaton\data\two_signals\time_evolution\III_chaotic_scan_p_ini';

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
p1 = 0:0.1:1;
p2 = 0:0.1:1;
nruns = 10;
tmax = 10000;
filecount = zeros(numel(p1), numel(p2));
t_out_all = zeros(numel(p1), numel(p2), nruns); % final times
period = zeros(numel(p1), numel(p2), nruns); % periodicity test
t_onset = tmax*ones(numel(p1), numel(p2), nruns); % default set to tmax

for i=1:numel(names)
    disp(names{i});
    load(fullfile(path, names{i}));
    p_ini = save_consts_struct.p_ini;
    idx1 = find(p_ini(1) == p1, 1);
    idx2 = find(p_ini(2) == p2, 1);
    filecount(idx1, idx2) = filecount(idx1, idx2) + 1;
    idx3 = filecount(idx1, idx2);
    t_out_all(idx1, idx2, idx3) = t_out;
    [period(idx1, idx2, idx3), t_onset(idx1, idx2, idx3)] = ...
        periodicity_test_short(cells_hist, save_consts_struct); % periodicity test
end

%% Save the loaded data
fname_str = 'III_chaotic_scan_p_ini_data';
save_path = 'H:\My Documents\Multicellular automaton\data\two_signals\time_evolution';
%save_vars = {'p1', 'p2', 'filecount', 't_out_all', 'period'};
save(fullfile(save_path, strcat(fname_str, '.mat') ), 'p1', 'p2', 'filecount', 't_out_all', 'period') %, save_vars);
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
xlabel('$$p_1$$')
ylabel('$$p_2$$')
set(gca, 'FontSize', 20);
ylim(c, [0 1]);
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
%ylim(c, [0 1]);

%% Analyze periods
% distribution of periods
h4 = figure(4);
histogram(period(period~=0) );
xlabel('Period');
ylabel('Count');
set(gca, 'FontSize', 20);

%% average period (count only found periods)
period_mean = zeros(numel(p1), numel(p2));
period_count = zeros(numel(p1), numel(p2));
for i=1:numel(t_out_all)
    if period(i)>0
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
xlabel('$$p_1$$')
ylabel('$$p_2$$')
ylabel(c, '$$\langle \tau \rangle$$', ...
    'Interpreter', 'latex', 'FontSize', 20);
set(gca, 'FontSize', 20);