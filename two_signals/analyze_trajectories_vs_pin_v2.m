%% Analyze saved trajectories across a range of pin
% v2: for trajectories that already have calculated periods (short style)
clear all
close all
set(0, 'defaulttextinterpreter', 'latex');
%%
%{
path = 'H:\My Documents\Multicellular automaton\data\two_signals\time_evolution\vs_p_ini_batch3';
%path = 'D:\Multicellularity\data\two_signals\time_evolution\vs_pini_batch3';
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

%% Load raw data
N = 225;
nruns = 80;
tmax = 10000;

p1 = 0:0.1:1;
p2 = 0:0.1:1;
N1 = round(p1*N);
N2 = round(p2*N);

filecount = zeros(numel(p1), numel(p2));
t_out_all = zeros(numel(p1), numel(p2), nruns); % final times
period_all = zeros(numel(p1), numel(p2), nruns); % periodicity test
t_onset_all = zeros(numel(p1), numel(p2), nruns); 

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
            period_all(idx1, idx2, idx3) = period;
            t_onset_all(idx1, idx2, idx3) = t_onset;
        end
    end
end
%}
%% Save the loaded raw data
%{
%fname_str = sprintf('III_chaotic_scan_p_ini_batch2_data_nruns%d', nruns);
fname_str = sprintf('vs_p_ini_batch3_data_nruns%d', nruns);
save_path = 'H:\My Documents\Multicellular automaton\data\two_signals\time_evolution';
%save_path = 'D:\Multicellularity\data\two_signals\time_evolution';

%save_vars = {'p1', 'p2', 'filecount', 't_out_all', 'period'};
%save(fullfile(save_path, strcat(fname_str, '.mat') ), 'p1', 'p2', 'filecount', 't_out_all', 'period', 't_onset') %, save_vars);
save(fullfile(save_path, strcat(fname_str, '.mat') ), 'p1', 'p2', 'filecount', 't_out_all', 'period_all', 't_onset_all') %, save_vars);

% folder for saving figures
save_path_fig = 'H:\My Documents\Multicellular automaton\figures\two_signals\analyze_trajectories_vs_pin';
%}
%% Load processed data
N = 225;
nruns = 80;
tmax = 10000;
%fname_str = 'III_chaotic_scan_p_ini_data';
fname_str = sprintf('vs_p_ini_batch3_data_nruns%d', nruns);
load_path = 'H:\My Documents\Multicellular automaton\data\two_signals\time_evolution';
load(fullfile(load_path, strcat(fname_str, '.mat')));

% folder for saving figures
save_path_fig = 'H:\My Documents\Multicellular automaton\figures\two_signals\analyze_trajectories_vs_pin';

%% Processing
% set period of non-periodic trajectories to 0
period_all(((period_all==Inf) + (t_out_all<tmax))==2) = 0;
%for i=1:numel(period_all)
%    if period_all(i)==Inf && t_out_all(i)<tmax
%        period_all(i) = 0;
%    end
%end

class_all = zeros(size(period_all));
for i=1:numel(class_all)
    if period_all(i)==0 % non-period, regular
        class_all(i) = 1;
    elseif period_all(i)==Inf % chaotic
        class_all(i) = 3;
    else % periodic
        class_all(i) = 2;
    end
end

%% Classify trajectories
C = categorical(class_all, 1:3, {'regular', 'periodic', 'chaotic'});

% Pie chart
h1=figure(1);
p = pie(C);
p(2).FontSize = 24;
p(4).FontSize = 24;
title(sprintf('%d simulations', numel(class_all)), 'FontSize', 24);

% modify text labels
ax = gca;
ax.Children(1).String = sprintf('periodic: %d', sum(class_all(:)==2));
ax.Children(3).String = sprintf('regular: %d', sum(class_all(:)==1));

qsave = 1;
if qsave
    fname = fullfile(save_path_fig, strcat(fname_str, '_class_pie_chart'));
    save_figure(h1, 10, 8, fname, '.pdf');
end
%% Classification vs p_ini
%frac_regular = zeros(numel(p1), numel(p2));
frac_regular  = sum(class_all == 1, 3)/nruns;
frac_periodic  = sum(class_all == 2, 3)/nruns;
frac_chaotic  = sum(class_all == 3, 3)/nruns;

% Plot fraction "regular"
h21=figure(21);
imagesc(p1, p2, frac_regular);
set(gca, 'YDir', 'Normal');
c = colorbar;
xlabel('$$p_1$$')
ylabel('$$p_2$$')
ylabel(c, 'Fraction regular', ...
    'Interpreter', 'latex', 'FontSize', 20);
set(gca, 'FontSize', 20);
caxis([0 1]);
title('Regular');
%ylim(c, [0 tmax]);

qsave = 1;
if qsave
    fname = fullfile(save_path_fig, strcat(fname_str, '_frac_regular_vs_p1_p2'));
    save_figure(h21, 10, 8, fname, '.pdf');
end

% Plot fraction periodic
h22=figure(22);
imagesc(p1, p2, frac_periodic);
set(gca, 'YDir', 'Normal');
c = colorbar;
xlabel('$$p_1$$')
ylabel('$$p_2$$')
ylabel(c, 'Fraction periodic', ...
    'Interpreter', 'latex', 'FontSize', 20);
set(gca, 'FontSize', 20);
caxis([0 1]);
title('Periodic');
%ylim(c, [0 tmax]);

qsave = 1;
if qsave
    fname = fullfile(save_path_fig, strcat(fname_str, '_frac_periodic_vs_p1_p2'));
    save_figure(h22, 10, 8, fname, '.pdf');
end

% Plot fraction chaotic
h23=figure(23);
imagesc(p1, p2, frac_chaotic);
set(gca, 'YDir', 'Normal');
c = colorbar;
xlabel('$$p_1$$')
ylabel('$$p_2$$')
ylabel(c, 'Fraction chaotic', ...
    'Interpreter', 'latex', 'FontSize', 20);
set(gca, 'FontSize', 20);
caxis([0 1]);
title('Chaotic');
%ylim(c, [0 tmax]);

qsave = 1;
if qsave
    fname = fullfile(save_path_fig, strcat(fname_str, '_frac_chaotic_vs_p1_p2'));
    save_figure(h23, 10, 8, fname, '.pdf');
end

%% Distribution of final t
h3=figure(3);
bins = 0:tmax/30:tmax;
histogram(t_out_all, bins);
xlabel('Final time');
ylabel('Count');
%title(sprintf('%d trajectories, %d periodic, %d chaotic',...
%    sum(sum(filecount)), n_periodic, n_chaotic));
title('All trajectories');
set(gca, 'FontSize', 20);

qsave = 1;
if qsave
    fname = fullfile(save_path_fig, strcat(fname_str, '_t_out_hist'));
    save_figure(h3, 10, 8, fname, '.pdf');
end

%%
% over regular trajectories
h31=figure(31);
bins = 0:tmax/30:tmax;
histogram(t_out_all(class_all==1), bins);
xlabel('Final time');
ylabel('Count');
title('Regular');
set(gca, 'FontSize', 20);

qsave = 1;
if qsave
    fname = fullfile(save_path_fig, strcat(fname_str, '_t_out_regular_hist'));
    save_figure(h31, 10, 8, fname, '.pdf');
end

% over periodic trajectories
h32=figure(32);
bins = 0:tmax/30:tmax;
histogram(t_onset_all(class_all==2), bins);
xlabel('Time onset periodicity');
ylabel('Count');
title('Periodic');
set(gca, 'FontSize', 20);

qsave = 1;
if qsave
    fname = fullfile(save_path_fig, strcat(fname_str, '_t_out_periodic_hist'));
    save_figure(h32, 10, 8, fname, '.pdf');
end

%% t_out vs p_ini
% Calculations
t_out_mean_all = zeros(numel(p1), numel(p2));
t_out_mean_regular = zeros(numel(p1), numel(p2)); %regular
%regular_count = t_out_mean_regular; % count # regular
t_out_mean_periodic = zeros(numel(p1), numel(p2)); %periodic
%periodic_count = t_out_mean_periodic; % count # regular

%t_out_all((period_all==0)+(t_out_all<tmax)==2);
%sum((period_all==0)+(t_out_all<tmax)==1, 3);

for i=1:numel(t_out_all)
    [i1,i2,i3] = ind2sub(size(t_out_all), i); 
    this_t_out = t_out_all(i1,i2,i3);
    t_out_mean_all(i1, i2) = t_out_mean_all(i1, i2) + this_t_out;
    if t_out_all(i)<tmax 
        if period_all(i)==0 % regular
            t_out_mean_regular(i1, i2) = t_out_mean_regular(i1, i2) + this_t_out;
            %regular_count(i1, i2) = regular_count(i1, i2) + 1;
        else % periodic
            t_out_mean_periodic(i1, i2) = t_out_mean_periodic(i1, i2) + this_t_out;
            %periodic_count(i1, i2) = periodic_count(i1, i2) + 1;
        end
    end
end

% recall:
%frac_regular  = sum(class_all == 1, 3)/nruns;
%frac_periodic  = sum(class_all == 2, 3)/nruns;
t_out_mean_all = t_out_mean_all./nruns;
t_out_mean_regular = t_out_mean_regular./(frac_regular*nruns);
t_out_mean_periodic = t_out_mean_periodic./(frac_periodic*nruns);
%}
%%
% <t_out> over all trajectories
h4 = figure(4);
imagesc(p1, p2, t_out_mean_all);
set(gca, 'YDir', 'Normal');
c = colorbar;
xlabel('$$p_1$$')
ylabel('$$p_2$$')
ylabel(c, '$$\langle t_{out} \rangle$$', ...
    'Interpreter', 'latex', 'FontSize', 20);
set(gca, 'FontSize', 20);
dc = 500; % round up to nearest multiple of dc
cmax = ceil(max(t_out_mean_all(:))/dc)*dc;
caxis([0 cmax]);
title('All trajectories');

qsave = 1;
if qsave
    fname = fullfile(save_path_fig, strcat(fname_str, '_mean_t_out_vs_p1_p2'));
    save_figure(h4, 10, 8, fname, '.pdf');
end
%% <t_out> over regular trajectories
h41 = figure(41);
imagesc(p1, p2, t_out_mean_regular);
set(gca, 'YDir', 'Normal');
c = colorbar;
xlabel('$$p_1$$')
ylabel('$$p_2$$')
ylabel(c, '$$\langle t_{out} \rangle$$ (regular)', ...
    'Interpreter', 'latex', 'FontSize', 20);
set(gca, 'FontSize', 20);
title('Regular');
%caxis([0 tmax]);
%ylim(c, [0 tmax]);
dc = 200; % round up to nearest multiple of dc
cmax = ceil(max(t_out_mean_regular(:))/dc)*dc;
caxis([0 cmax]);
%ylim(c, [0 tmax]);

qsave = 1;
if qsave
    fname = fullfile(save_path_fig, strcat(fname_str, '_mean_t_out_vs_p1_p2_regular'));
    save_figure(h41, 10, 8, fname, '.pdf');
end
%%
% <t_out> over periodic trajectories
h42 = figure(42);
imagesc(p1, p2, t_out_mean_periodic);
set(gca, 'YDir', 'Normal');
c = colorbar;
xlabel('$$p_1$$')
ylabel('$$p_2$$')
ylabel(c, '$$\langle t_{out} \rangle$$ (periodic)', ...
    'Interpreter', 'latex', 'FontSize', 20);
set(gca, 'FontSize', 20);
title('Periodic');
dc = 500; % round up to nearest multiple of dc
cmax = ceil(max(t_out_mean_periodic(:))/dc)*dc;
caxis([0 cmax]);
%ylim(c, [0 tmax]);

qsave = 1;
if qsave
    fname = fullfile(save_path_fig, strcat(fname_str, '_mean_t_out_vs_p1_p2_periodic'));
    save_figure(h42, 10, 8, fname, '.pdf');
end

%% Analyze periods
% number of chaotic trajectories
n_chaotic = sum((period_all(:)==Inf).*(t_out_all(:) == tmax));

% data
period_data = period_all(period_all~=0);
uniq = unique(period_data);
C = categorical(period_data, uniq);

n_periodic = sum(period_data(:)~=Inf);

% distribution of periods
h5 = figure(5);
histogram(C);
xlabel('Period');
ylabel('Count');
title(sprintf('%d trajectories, %d periodic, %d chaotic',...
    sum(sum(filecount)), n_periodic, n_chaotic));
set(gca, 'FontSize', 20);

qsave = 1;
if qsave
    fname = fullfile(save_path_fig, strcat(fname_str, '_period_hist_v2'));
    save_figure(h5, 12, 8, fname, '.pdf');
end

%% Average period (periodic trajectories)
period_mean = zeros(numel(p1), numel(p2));
period_count = zeros(numel(p1), numel(p2));
for i=1:numel(t_out_all)
    if period_all(i)>0 && period_all(i)<Inf
        [i1,i2,i3] = ind2sub(size(period_all), i); 
        period_count(i1, i2) = period_count(i1, i2) + 1;
        period_mean(i1, i2) = period_mean(i1, i2) + period_all(i1,i2,i3);
    end
end
period_mean = period_mean./period_count;

h6 = figure(6);
imagesc(p1, p2, period_mean);
set(gca, 'YDir', 'Normal');
c = colorbar;
caxis([0 max(uniq(uniq<Inf))])
xlabel('$$p_1$$')
ylabel('$$p_2$$')
ylabel(c, 'Mean period (time steps)', ...
    'Interpreter', 'latex', 'FontSize', 20);
set(gca, 'FontSize', 20);
dc = 5; % round up to nearest multiple of dc
cmax = ceil(max(period_mean(:))/dc)*dc;
caxis([0 cmax]);

qsave = 1;
if qsave
    fname = fullfile(save_path_fig, strcat(fname_str, '_mean_period_vs_p1_p2'));
    save_figure(h6, 10, 8, fname, '.pdf');
end