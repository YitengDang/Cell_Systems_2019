%% Analyze saved trajectories across a range of pin
% v3: for analyzing travelling waves
clear all
close all
set(0, 'defaulttextinterpreter', 'latex');
%% Parameters
N = 225;
a0 = 1.5;
nruns = 100;
tmax = 10000;

p1_all = 0:0.1:1;
p2_all = 0:0.1:1;
N1 = round(p1_all*N);
N2 = round(p2_all*N);

K12 = 22;

% folder for saving figures
save_path_fig = 'H:\My Documents\Multicellular automaton\figures\two_signals\trajectories_vs_pin\trav_wave';
%% Get all filenames
%
parent_folder = fullfile('N:\tnw\BN\HY\Shared\Yiteng\two_signals\travelling_wave_analysis',...
    sprintf('vs_p0_K12_%d', K12));

names = cell(numel(p1), numel(p2), nruns);
for i1=1:numel(p1)
    for i2=1:numel(p2)
        p = [p1(i1) p2(i2)];
        subfolder = strrep(sprintf('ini_p1_%.2f_p2_%.2f', p(1), p(2)), '.', 'p');
        folder = fullfile(parent_folder, subfolder);
        
        listing = dir(folder);
        num_files = numel(listing)-2; %first two entries are not useful
        count = 0;
        for i = 1:num_files
            filename = listing(i+2).name;
            % remove extension and do not include txt files
            [~,fname,ext] = fileparts(filename);
            if strcmp(ext, '.mat')
                count = count + 1;
                names{i1, i2, count} = fname;
            end
        end
    end
end

%% Load raw data
% possible filter to use
pattern = '.'; % '.' = anything 'p_ini\dp\d{2}_0p30';

% variables to store
filecount = zeros(numel(p1), numel(p2));
t_out_all = zeros(numel(p1), numel(p2), nruns); % final times
period_all = zeros(numel(p1), numel(p2), nruns); % periodicity test
t_onset_all = zeros(numel(p1), numel(p2), nruns); 
trav_wave_all = zeros(numel(p1), numel(p2), nruns); 
trav_wave_2_all = zeros(numel(p1), numel(p2), nruns); 

error_files = {};
for i1=1:numel(p1)
    for i2=1:numel(p2)
        p = [p1(i1) p2(i2)];
        subfolder = strrep(sprintf('ini_p1_%.2f_p2_%.2f', p(1), p(2)), '.', 'p');
        folder = fullfile(parent_folder, subfolder);
        for i3=1:nruns
            fname = names{i1, i2, i3};
            if isempty(regexp(fname, pattern, 'once')) % only load files matching a certain pattern
                continue
            else
                disp(fname);
                load(fullfile(folder, fname ));
                
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
                    
                    trav_wave_all(idx1, idx2, idx3) = trav_wave;
                    
                    [trav_wave_new, trav_wave_2] = travelling_wave_test(cells_hist, a0, period, t_out);
                    if trav_wave_new~=trav_wave
                        warning('Travelling wave tests do not match!');
                        error_files{end+1} = fname;
                    end
                    trav_wave_2_all(idx1, idx2, idx3) = trav_wave_2;
                end
            end
        end
    end
end
%}
%% Save the loaded raw data
%
fname_str = sprintf('trav_wave_occur_vs_p_ini_set2_K12_9_nruns%d', nruns);
save_path = fullfile('N:\tnw\BN\HY\Shared\Yiteng\two_signals\travelling_wave_analysis',...
    sprintf('vs_p0_K12_%d', K12));

save(fullfile(save_path, strcat(fname_str, '.mat') ), 'p1', 'p2', 'filecount',...
    't_out_all', 'period_all', 't_onset_all', 'trav_wave_all', 'trav_wave_2_all') %, save_vars);
%}
%% Load processed data
fname_str = sprintf('trav_wave_occur_vs_p_ini_set2_K12_%d_nruns%d', K12, nruns);
load_path = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\travelling_wave_analysis';
subfolder = sprintf('vs_p0_K12_%d', K12);
load( fullfile(load_path, subfolder, strcat(fname_str, '.mat')) );

%% Processing: classify trajectories
% class 1: non-periodic
% class 2: periodic, not travelling wave
% class 3: periodic travelling wave
% class 4: tmax reached, inconclusive

% store all classes
class_all = zeros(size(period_all));
class_all_2 = zeros(size(period_all)); % based on less stringent trav. wave criterion
n_classes = 4;

% set period of non-periodic trajectories to 0
idx_class = {};
idx_class{1} = ((period_all==Inf) + (t_out_all<tmax))==2;
idx_class{2} = ((period_all<Inf) + (trav_wave_all==0))==2;
idx_class{3} = ((period_all<Inf) + (trav_wave_all==1))==2;
idx_class{4} = ((period_all==Inf) + (t_out_all==tmax))==2;

idx_class_2 = {}; 
idx_class_2{2} = ((period_all<Inf) + (trav_wave_2_all==0))==2;
idx_class_2{3} = ((period_all<Inf) + (trav_wave_2_all==1))==2;

class_all(idx_class{1}) = 1;
class_all(idx_class{2}) = 2;
class_all(idx_class{3}) = 3;
class_all(idx_class{4}) = 4;

class_all_2(idx_class{1}) = 1;
class_all_2(idx_class_2{2}) = 2;
class_all_2(idx_class_2{3}) = 3;
class_all_2(idx_class{4}) = 4;

% Classify trajectories
%C = categorical(class_all, 1:n_classes); %, {'non-periodic', 'periodic, non-trav. wave', 'trav. wave', 'unknown'});
X = [sum(idx_class{1}(:)) sum(idx_class{2}(:)) sum(idx_class{3}(:)) sum(idx_class{4}(:))]; 
X2 = [sum(idx_class{1}(:)) sum(idx_class_2{2}(:)) sum(idx_class_2{3}(:)) sum(idx_class{4}(:))]; 

% Fractions
%frac_regular = zeros(numel(p1), numel(p2));
frac_regular = sum(class_all == 1, 3)/nruns;
frac_periodic = sum(((class_all == 2) + (class_all == 3)), 3)/nruns;
frac_trav_wave = sum(class_all == 3, 3)/nruns;
frac_trav_wave_2 = sum(class_all_2 == 3, 3)/nruns;

%% calculated weighted estimate of total fraction of each class
% N.B. weighted average is not over entire distribution (p1, p2 values), but
% over only a part of it. This shouldn't affect the results.
density_states = zeros(numel(p1_all), numel(p2_all));
sum1 = 0;
for i=1:numel(p1_all)
    for j=1:numel(p2_all)
        n1 = round(p1_all(i)*N);
        n2 = round(p2_all(j)*N);
        density_states(i,j) = (nchoosek(N, n1)/2^N)*(nchoosek(N, n2)/2^N);
    end
end
density_states = density_states/(sum(sum(density_states)));

X_count_by_p = zeros(4, numel(p1_all), numel(p2_all));
X_weighted = zeros(1, 4);
for i=1:4
    class_count = sum(idx_class{i}, 3);
    X_count_by_p(i,:,:) = class_count;
    X_weighted(i) = sum(sum( class_count/nruns.*density_states ));
end

%% Pie chart
h1=figure(1);
p = pie(X2);
%delete(p(2));
%delete(p(4));
%delete(p(6));
p(2).FontSize = 16;
p(4).FontSize = 16;
p(6).FontSize = 16;
%p(8).FontSize = 16;
title(sprintf('%d simulations', numel(class_all)), 'FontSize', 24);
legend({'non-periodic', 'periodic, non-trav. wave', 'trav. wave', 'unknown'},...
    'Location', 'so', 'Orientation','horizontal', 'FontSize', 20);

% modify text labels
%ax = gca;
%ax.Children(1).String = sprintf('periodic: %d', sum(class_all(:)==2));
%ax.Children(3).String = sprintf('regular: %d', sum(class_all(:)==1));

qsave = 0;
if qsave
    fname = fullfile(save_path_fig, strcat(fname_str, '_class_pie_chart_class_2'));
    save_figure(h1, 10, 8, fname, '.pdf');
end

%% Stacked bar chart for fractions of each type (regular, periodic, TW, ...)
% class 1: non-periodic
% class 2: periodic, not travelling wave
% class 3: periodic travelling wave
% class 4: tmax reached, inconclusive
% calculate fractions of simulations belonging to each class
class_frac = zeros(2, 4);
class_frac(1,:) = X/sum(X); % unweighted / raw
class_frac(2,:) = X_weighted;

h=figure;
bar(class_frac, 'stacked');
xlabel('');
ylabel('Fraction');
legend({'non-periodic', 'periodic, non-TW', 'travelling wave', 'inconclusive'}, ...
    'Location', 'eo');
set(gca, 'XTick', 1:2, 'XTickLabels', {'Unweighted', 'Weighted'})
set(gca, 'FontSize', 20);
%xlim([0 4]);

qsave = 1;
fname = fullfile(save_path_fig, strcat(fname_str, '_class_frac_overall'));
save_figure(h, 10, 8, fname, '.pdf', qsave);

%}
%% Classification vs p_ini
% Plot fraction "regular"
h21=figure(21);
imagesc(p1, p2, frac_regular);
set(gca, 'YDir', 'Normal');
c = colorbar;
xlabel('$$p_1$$')
ylabel('$$p_2$$')
ylabel(c, 'Fraction', ...
    'Interpreter', 'latex', 'FontSize', 20);
set(gca, 'FontSize', 20);
caxis([0 1]);
title('Non-periodic');
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
ylabel(c, 'Fraction', ...
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

% Plot fraction travelling waves
h23=figure(23);
imagesc(p1, p2, frac_trav_wave);
set(gca, 'YDir', 'Normal');
c = colorbar;
xlabel('$$p_1$$')
ylabel('$$p_2$$')
ylabel(c, 'Fraction', ...
    'Interpreter', 'latex', 'FontSize', 20);
set(gca, 'FontSize', 20);
caxis([0 1]);
title('Fraction trav. waves');
%ylim(c, [0 tmax]);

qsave = 1;
if qsave
    fname = fullfile(save_path_fig, strcat(fname_str, '_frac_trav_waves_vs_p1_p2'));
    save_figure(h23, 10, 8, fname, '.pdf');
end

%% Plot fraction travelling waves, class 2
h24=figure(24);
imagesc(p1, p2, frac_trav_wave_2);
set(gca, 'YDir', 'Normal');
c = colorbar;
xlabel('$$p_1$$')
ylabel('$$p_2$$')
ylabel(c, 'Fraction', ...
    'Interpreter', 'latex', 'FontSize', 20);
set(gca, 'FontSize', 20);
caxis([0 1]);
title('Fraction trav. waves');
%ylim(c, [0 tmax]);

qsave = 1;
if qsave
    fname = fullfile(save_path_fig, strcat(fname_str, '_frac_trav_waves_vs_p1_p2_class2'));
    save_figure(h24, 10, 8, fname, '.pdf');
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
histogram(t_onset_all(class_all==2 + class_all==3), bins);
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