%% Analyze saved trajectories for signs of travelling waves 
% Check across a range of a variable (= noise, Hill)
clear all
close all
set(0, 'defaulttextinterpreter', 'tex');
%% Parameters
N_all = [8 10 12 14 16 18 20 25].^2;
a0 = 1.5;
nruns = 300;
tmax = 10000;
%noise_all = [0 0.01 0.05 0.1 0.5];
%hill_all = [2 5 10 100 Inf];
K12 = 10;
digits = 5;

% variable to loop over
loopvar = N_all; %hill_all; %noise_all; %
%loopvar_str = 'hill'; %'noise'; 
loopvar_str = 'N';

%% Get all filenames
%
%path = 'H:\My Documents\Multicellular automaton\data\two_signals\time_evolution\vs_p_ini_batch3';
%path = 'D:\Multicellularity\data\two_signals\time_evolution\vs_pini_batch3';
%parent_folder = 'L:\BN\HY\Shared\Yiteng\two_signals\travelling_wave_analysis\vs_noise\K12_9';
%parent_folder = 'L:\BN\HY\Shared\Yiteng\two_signals\travelling_wave_analysis\vs_Hill\K12_9';
parent_folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\trav_wave_vs_N';

names = cell(numel(loopvar), nruns);
for i1=1:numel(loopvar)
    %folder = parent_folder; %fullfile(parent_folder, subfolder);
    subfolder = sprintf('N%d', N_all(i1));
    folder = fullfile(parent_folder, subfolder);
    listing = dir(folder);
    num_files = numel(listing)-2; %first two entries are not useful
    count = 0;
    for i2 = 1:nruns %num_files
        filename = listing(i2+2).name;
        % remove extension and do not include txt files
        [~,fname,ext] = fileparts(filename);
        if strcmp(ext, '.mat')
            count = count + 1;
            names{i1, count} = fname;
        end
    end
end
%}
%% filename pattern
% get file info from regexp
pat1 = '\d+|Inf';
pat2 = '(\d+p\d+|Inf)';
pat3 = '\w*';

if strcmp(loopvar_str, 'noise')
    pattern = sprintf('two_signal_mult_N%d_noise_%s_K12_%d_t_out_%s_period_%s-v%s%s',...
        N, pat2, K12, pat1, pat1, pat1, pat3);
    %pattern = sprintf('two_signal_mult_N%d_initiateI0_noise_%s_t_out_%s_period_%s_%stemp-v%s',...
    %    N, pat2, pat1, pat1, pat3, pat1); % noise
elseif strcmp(loopvar_str, 'hill')
    pattern = sprintf('two_signal_mult_N%d_hill_%s_K12_%d_t_out_%s_period_%s-v%s%s',...
        N, pat2, K12, pat1, pat1, pat1, pat3); % Hill
    %pattern = sprintf('two_signal_mult_N%d_initiateI0_hill_%s_t_out_%s_period_%s_%stemp-v%s',...
    %    N, pat2, pat1, pat1, pat3, pat1); % Hill
elseif strcmp(loopvar_str, 'N')
    pattern = sprintf('two_signal_mult_N%s_t_out_%s_period_%s-v%s',...
        '(\d+)', '\d+', '(\d+|Inf)', '\d+'); % Hill
    pattern2 = sprintf('two_signal_mult_N%s_initiateI0_randpos_mcsteps0_K12_10_t_out_%s_period_%s-v%s',...
        '(\d+)', '\d+', '(\d+|Inf)', '\d+');
end
%% Load raw data
%{
% variables to store
fname_all = cell(numel(loopvar), nruns);
filecount = zeros(numel(loopvar), 1);
t_out_all = zeros(numel(loopvar), nruns); % final times
period_all = zeros(numel(loopvar), nruns); % periodicity test
t_onset_all = zeros(numel(loopvar), nruns);
trav_wave_all = zeros(numel(loopvar), nruns);
trav_wave_2_all = zeros(numel(loopvar), nruns);

error_files = {};

for i1=1:numel(loopvar)
    for i2=1:nruns
        %----------------------------------------------------------------------
        %test_fname1 = 'two_signal_mult_N225_initiateI0_noise_0p50_t_out_10000_period_Inf_tmax_reached_temp-v1';
        %test_fname2 = 'two_signal_mult_N225_initiateI0_noise_0p10_t_out_3491_period_15_temp-v1';
        fname = names{i1, i2};

        [~, tokens] = regexp(fname, pattern, 'match', 'tokens');
        if isempty(tokens) % temporary fix
            [~, tokens] = regexp(fname, pattern2, 'match', 'tokens');
        end
        % get variable from regexp tokens, match to variable list
        this_var = str2double(strrep(tokens{1}{1}, 'p', '.')); % noise / Hill

        idx = find(this_var == loopvar);

        if isempty(idx) % only load files matching a certain pattern
            continue
        else
            disp(fname);
            
            subfolder = sprintf('N%d', N_all(idx));
            folder = fullfile(parent_folder, subfolder);

            load(fullfile(folder, fname));

            %{
            filecount(idx) = filecount(idx) + 1;
            idx2 = filecount(idx);
            fname_all{idx, idx2} = fname; % store filename
            t_out_all(idx, idx2) = t_out;
            period_all(idx, idx2) = period;
            t_onset_all(idx, idx2) = t_onset;
            %trav_wave_all(idx, idx2) = trav_wave;
            if period<Inf
                [trav_wave, trav_wave_2] = travelling_wave_test(cells_hist, a0,...
                    period, t_out, distances, digits);
                trav_wave_all(idx, idx2) = trav_wave;
                trav_wave_2_all(idx, idx2) = trav_wave_2;
            % else: not a travelling wave, no need to update trav_wave_all
            end
            %}
        end
    end
end
%}

%% tester
%{
folder = 'L:\BN\HY\Shared\Yiteng\two_signals\travelling_wave_analysis\vs_Hill';
fname = 'two_signal_mult_N225_initiateI0_hill_100p00_t_out_1503_period_15_temp-v1';
load(fullfile(folder, fname));
if period<Inf
    [trav_wave, trav_wave_2] = travelling_wave_test(cells_hist, a0, period, t_out);
end
%}
%{
test_str = fname; % 'two_signal_mult_N225_initiateI0_hill_5p00_t_out_9495_period_10_tmax_reached_temp-v1';
[match, tokens] = regexp(test_str, pattern, 'match', 'tokens');
this_var = str2double(strrep(tokens{1}{1}, 'p', '.')); % noise / Hill
%}
% Select data to examine
%{
%select_idx = (trav_wave_all ~= trav_wave_2_all);
select_idx = trav_wave_all==1;
fname_all(select_idx);
%}
%% Save the loaded raw data
%{
%fname_str = sprintf('trav_wave_occur_vs_%s_K12_%d_nruns%d_round_p_I_digits_%d',...
%    loopvar_str, K12, nruns, digits);
fname_str = sprintf('trav_wave_occur_vs_%s_K12_%d_nruns%d_round_p_I_digits_%d',...
    loopvar_str, K12, nruns, digits);
%save_path = 'L:\BN\HY\Shared\Yiteng\two_signals\travelling_wave_analysis\analyzed_data';
save_path = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\trav_wave_vs_N';

save(fullfile(save_path, strcat(fname_str, '.mat') ), 'loopvar', 'loopvar_str', 'filecount',...
    't_out_all', 'period_all', 't_onset_all', 'trav_wave_all', 'trav_wave_2_all') %, save_vars);
%}

%% Load processed data
%
% vs noise, Hill
%fname_str = sprintf('trav_wave_occur_vs_%s_K12_%d_nruns%d_round_p_I_digits_%d',...
%    loopvar_str, K12, nruns, digits);
%load_path = 'L:\BN\HY\Shared\Yiteng\two_signals\travelling_wave_analysis\analyzed_data';
% vs N
load_path = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\trav_wave_vs_N';
fname_str = 'trav_wave_occur_vs_N_K12_10_nruns300_round_p_I_digits_5';

load(fullfile(load_path, strcat(fname_str, '.mat')));
%}

%% Fraction travelling waves
% folder for saving figures
save_path_fig = 'H:\My Documents\Multicellular automaton\figures\two_signals\trav_wave_vs_N';

h = figure;
hold on
%---
if strcmp(loopvar_str, 'hill')
    xdata = [loopvar(1:end-1) 10^3];
    plot(xdata, sum(trav_wave_all, 2)/nruns, 'bo--', 'LineWidth', 2);
    plot(xdata, sum(trav_wave_2_all, 2)/nruns, 'ro--', 'LineWidth', 2);
    xlabel('Hill coeff. $n$')
    xtick_str = [string(loopvar(1:end-1)), "Inf"];
    set(gca, 'XScale', 'log');
%---
elseif strcmp(loopvar_str, 'noise')
    xdata = [0.001 loopvar(2:end)];
    plot(xdata, sum(trav_wave_all, 2)/nruns, 'bo--', 'LineWidth', 2);
    plot(xdata, sum(trav_wave_2_all, 2)/nruns, 'ro--', 'LineWidth', 2);
    xlabel('Noise $\alpha$')
    xtick_str = string([0 loopvar(2:end)]);
    set(gca, 'XScale', 'log');
%---
elseif strcmp(loopvar_str, 'N')
    xdata = loopvar;
    plot(xdata, sum(trav_wave_all, 2)/nruns, 'bo--', 'LineWidth', 2);
    plot(xdata, sum(trav_wave_2_all, 2)/nruns, 'ro--', 'LineWidth', 2);
    xlabel('$N$')
    xtick_str = string(loopvar);
end
%---
ylabel('Fraction')
title('Fraction travelling waves');
set(gca, 'FontSize', 24);
set(gca, 'XTick', xdata, 'XTickLabels',...
    xtick_str);
ylim([0 1]);
legend('strict crit.', 'loose crit.');

qsave = 0;
if qsave
    fname = fullfile(save_path_fig, sprintf('frac_trav_waves_vs_%s_K12_%d_nruns_%d_p_I_%d_digits',...
        loopvar_str, K12, nruns, digits));
    save_figure(h, 10, 8, fname, '.pdf', qsave);
end

%% Filter data selectively
%{
% multiples of period 15
frac_period_15_mult = sum((mod(period_all, sqrt(N)) == 0), 2)/nruns;

% Estimate by hand
%frac_trav_wave_manual = [0 4 33 89 69]/100;

h = figure;
hold on
if strcmp(loopvar_str, 'hill')
    %---
    xdata = [loopvar(1:end-1) 10^3];
    xlabel('Hill coeff. $n$')
    xtick_str = [string(loopvar(1:end-1)), "Inf"];
    %---
elseif strcmp(loopvar_str, 'noise')
    xdata = [0.001 loopvar(2:end)];
    xlabel('Noise $\alpha$')
    xtick_str = string([0 loopvar(2:end)]);
end
%---
plot(xdata, frac_period_15_mult, 'go--', 'LineWidth', 2);
%plot(xdata, frac_trav_wave_manual, 'mo--', 'LineWidth', 2);
ylabel('Fraction')
title('Fraction travelling waves');
set(gca, 'FontSize', 24);
set(gca, 'XScale', 'log', 'XTick', xdata, 'XTickLabels',...
    xtick_str);
ylim([0 1]);
legend('period = multiple of 15', 'Location', 'nw');

qsave = 1;
if qsave
    fname = fullfile(save_path_fig, sprintf(...
        'frac_trav_waves_vs_%s_method2_K12_%d_nruns_%d_p_I_%d_digits',...
        loopvar_str, K12, nruns, digits));
    save_figure(h, 10, 8, fname, '.pdf', qsave);
end
%}
%% Plot t_onset vs loopvar
% process data
t_onset_mean = Inf*ones(numel(loopvar), 1);
t_onset_std = zeros(numel(loopvar), 1);
for i=1:numel(loopvar)
    %idx = find(mod(period_all(i,:), sqrt(N)) == 0);
    %idx = find(trav_wave_all(i, :));
    idx = find(trav_wave_2_all(i, :));
    if ~isempty(idx)
        t_onset_mean(i) = mean( t_onset_all(i, idx) );
        t_onset_std(i) = std( t_onset_all(i, idx) );
    end
end

% calculate confidence interval
% t_onset_error_margin = zeros(numel(loopvar), 1);
z = 1.96;
trav_wave_count = sum(trav_wave_all, 2);
t_onset_error_margin = z*t_onset_std./sqrt( trav_wave_count );

if strcmp(loopvar_str, 'hill')
    %------
    xdata = [loopvar(1:end-1) 10^3];
    xtick_str = [string(loopvar(1:end-1)), "Inf"];
    xlabel_text = 'Hill coeff. $n$';
    xlim([1 10^3]);
    set(gca, 'XScale', 'log');
    %------
elseif strcmp(loopvar_str, 'noise')
    xdata = [0.001 loopvar(2:end)];
    xtick_str = string([0 loopvar(2:end)]);
    xlabel_text = 'Noise $\alpha$';
    xlim([0.001 1]);
    set(gca, 'XScale', 'log');
    %------
elseif strcmp(loopvar_str, 'N')
    xdata = sqrt(loopvar);
    %xlabel('$N$')
    xlabel_text = 'Grid size';
    xtick_str = string(sqrt(loopvar));  
end

%% calculate fit function
[ffit, S] = polyfit(xdata, t_onset_mean', 1);
ffit_ydata = ffit(1)*xdata + ffit(2);

%mdl = fitlm(xdata, t_onset_mean, 'linear');

%% Plot variability vs variable

h = figure;
plot(xdata, t_onset_std./t_onset_mean, 'o',  'MarkerSize', 10)
ylabel('Coefficient of variation \sigma/\mu');
ylim([0 1.2]);
%%
h = figure;
hold on
%plot(xdata, t_onset_mean,...
%    'g^', 'LineWidth', 2, 'MarkerSize', 10);
errorbar(xdata, t_onset_mean, t_onset_error_margin,...
    'g^', 'LineWidth', 2, 'MarkerSize', 10);
plot(xdata, ffit_ydata, 'g--', 'LineWidth', 2);
box on
xlabel(xlabel_text);
ylabel('Travelling wave formation time');
title('Travelling wave formation time');
set(gca, 'FontSize', 24);
set(gca, 'XTick', xdata, 'XTickLabels',...
    xtick_str);
xlim([xdata(1)-1 xdata(end)+1]);
ylim([0 7000]);

% Plot CV
yyaxis right
plot(xdata, t_onset_std./t_onset_mean, 'o',  'MarkerSize', 10)
ylabel('Coefficient of variation \sigma/\mu');

%legend('period = multiple of 15', 'manual', 'Location', 'nw');

qsave = 0;
if qsave
    fname = fullfile(save_path_fig, sprintf(...
        't_onset_trav_waves_vs_%s_K12_%d_nruns_%d_digits_%d_errorbar',...
        loopvar_str, K12, nruns, digits));
    save_figure(h, 10, 8, fname, '.pdf', qsave);
end

%% Plot together with theoretical predictions
% Load data
folder = 'H:\My Documents\Multicellular automaton\data\two_signals';
fname_str = 'NwaveDensity'; %'TravWaveDensity';
fname = fullfile(folder, fname_str);
xls_data = xlsread(fname);

% Plot theoretical estimates
idx_sel = 1:size(xls_data, 2);
%h = figure;
plot(xls_data(1,idx_sel), 1./xls_data(2,idx_sel), 'bo',...
    'LineWidth', 2, 'MarkerSize', 10);
set(gca, 'YScale', 'log');
legend({'Model simulations', 'Random process'}, 'Location', 'nw');
ylim([1 10^150]);

qsave = 0;
if qsave
    fname = fullfile(save_path_fig, sprintf(...
        't_onset_trav_waves_vs_%s_K12_%d_nruns_%d_digits_%d_sim_and_calc_for_main_fig',...
        loopvar_str, K12, nruns, digits));
    save_figure(h, 10, 8, fname, '.pdf', qsave);
end
% Box plot
%{
idx = (mod(period_all, sqrt(N)) == 0);
box_data = t_onset_all(idx);
group_data = repmat(1:5', 1, nruns);
group_data = group_data(idx);
boxplot(box_data, group_data, 'Labels', [string(loopvar(2:end-1)), "Inf"]);
%}
%% Plot t_onset histograms - Normalized onset times (together)
bins = 10;
edges = 0:1:30;
N_hist = zeros(numel(loopvar), numel(edges)-1 );

for ii=1:numel(loopvar)
    idx = find(trav_wave_2_all(ii, :));
    if ~isempty(idx)
        N_hist(ii,:) = histcounts(t_onset_all(ii, idx)./t_onset_mean, edges,...
            'normalization', 'probability');
    end
end

figure;
hold on
bincenters = (edges(1:end-1)+edges(2:end))/2;
plot( bincenters, N_hist, 'x-');
%set(gca, 'YScale', 'log')
legend( sprintfc('N=%d', loopvar) );

% (2) plot separate histograms 
edges = 0:0.5:5;
N_hist = zeros(numel(loopvar), numel(edges)-1 );

for i=1:numel(loopvar)
    idx = find(trav_wave_2_all(i, :));
    if ~isempty(idx)
        N_hist(i, :) = histcounts( t_onset_all(i, idx)/t_onset_mean(i), edges );
    end
end

edgecenters = (edges(1:end-1) + edges(2:end))/2;
%figure;
%hold on
for ii=1:numel(loopvar)
    figure;
    bar( edgecenters, N_hist(ii,:) );
end

%% Plot t_onset histograms - Unnormalized onset times (together)
h = figure;
hold on
nbins = 10;
edges = 0:1000:10000;
N_hist = zeros(numel(loopvar), numel(edges)-1 );
for ii=1:numel(loopvar)
    idx = find(trav_wave_2_all(ii, :));
    if ~isempty(idx)
        %histogram(t_onset_all(ii, idx), nbins);
        N_hist(ii,:) = histcounts(t_onset_all(ii, idx), edges,...
            'normalization', 'probability');
    end
end

bincenters = (edges(1:end-1)+edges(2:end))/2;
%waterfall(bincenters, sqrt(loopvar), N_hist)
%ribbon( N_hist')
b=bar3( N_hist', 1 );
set(b,'FaceAlpha', 0.5)
set(gca, 'XTick', 1:size(N_hist, 1), 'XTickLabel', sqrt(loopvar));
set(gca, 'YTick', 0.5:2:size(N_hist, 2)+0.5, 'YTickLabel', edges(1:2:end)); % 1:size(N_hist, 2),
ylabel('TW formation time');
xlabel('Grid Size');
zlabel('Probability');
set(gca, 'FontSize', 20);
view(45, 45);

folder_out = 'H:\My Documents\Multicellular automaton\figures\two_signals\trav_wave_vs_N';
fname_str_out = sprintf('t_onset_TW_distribution_all_K12_%d_nruns_%d_digits_%d_v2_touching_bars', K12, nruns, digits);
fname = fullfile(folder_out, fname_str_out);
set(h, 'Units', 'inches', 'Position', [1 1 16 9]);

qsave = 1;
save_figure(h, 16, 9, fname, '.pdf', qsave);

%% Plot t_onset histograms - unnormalized, separate

loopvar_idx = 5;
N = loopvar(loopvar_idx);
fprintf('N = %d \n', N);

h=figure;
idx = find(trav_wave_2_all(loopvar_idx, :));
edges = 0:500:5000;
histogram(t_onset_all(loopvar_idx, idx), edges, 'Normalization', 'probability');
xlabel('TW formation time');
ylabel('Probability');
set(gca, 'FontSize', 20);

folder_out = 'H:\My Documents\Multicellular automaton\figures\two_signals\trav_wave_vs_N';
fname_str_out = sprintf('t_onset_TW_distribution_single_N_%d_K12_%d_nruns_%d_digits_%d',...
    N, K12, nruns, digits);
fname = fullfile(folder_out, fname_str_out);

qsave = 1;
save_figure(h, 10, 8, fname, '.pdf', qsave);


%% Plot t_onset full data
if strcmp(loopvar_str, 'hill')
    %------
    xdata = [loopvar(1:end-1) 10^3];
    xtick_str = [string(loopvar(1:end-1)), "Inf"];
    xlabel_text = 'Hill coeff. $n$';
    xlim([1 10^3]);
    set(gca, 'XScale', 'log');
    %------
elseif strcmp(loopvar_str, 'noise')
    xdata = [0.001 loopvar(2:end)];
    xtick_str = string([0 loopvar(2:end)]);
    xlabel_text = 'Noise $\alpha$';
    xlim([0.001 1]);
    set(gca, 'XScale', 'log');
    %------
elseif strcmp(loopvar_str, 'N')
    xdata = sqrt(loopvar);
    xlabel_text = 'Grid size ($\sqrt{N}$)';
    xtick_str = string(sqrt(loopvar));  
end

h = figure;
hold on
max_t_onset = zeros(numel(loopvar), 1);
for i=1:numel(loopvar)
    idx = find(trav_wave_all(i, :));
    if ~isempty(idx)
        scatter(xdata(i)*ones(size(idx)), t_onset_all(i, idx), 75 );
        max_t_onset(i) = max(t_onset_all(i, idx));
    end
end
xlabel(xlabel_text);
ylabel('$t_{onset}$');
title('Onset time of trav. wave');
set(gca, 'FontSize', 24);
set(gca, 'XTick', xdata, 'XTickLabels',...
    xtick_str);

% labels above
%
labels = sprintfc('#TW=%d', sum(trav_wave_all, 2) );
xt = get(gca, 'XTick');
text(xt, max_t_onset+30, labels, 'FontSize', 20,...
    'HorizontalAlignment', 'center', 'VerticalAlignment','bottom')
%}

qsave = 1;
fname = fullfile(save_path_fig, sprintf(...
    't_onset_trav_waves_vs_%s_K12_%d_nruns_%d_digits_%d_scatter',...
    loopvar_str, K12, nruns, digits));
save_figure(h, 10, 8, fname, '.pdf', qsave);
%% Plot as box plot
h = figure;
hold on
data = cell(numel(loopvar), 1);
data_array = [];
g = [];
for i=1:numel(loopvar)
    idx = find(trav_wave_all(i, :));
    disp(numel(idx));
    if ~isempty(idx)
        data{i} = t_onset_all(i, idx);
        data_array = [data_array, t_onset_all(i, idx)];
        g = [g; i*ones(numel(idx), 1)];
    end
end
boxplot([data{:}], g, 'Labels', sprintfc('%d', sqrt(N_all)))
xlabel('Grid size ($\sqrt{N}$)');
%xlabel('$N$');
ylabel('$t_{onset}$');
set(gca, 'FontSize', 24);
% labels above
%
%labels = sprintfc('#TW=%d', sum(trav_wave_all, 2) );
labels = sprintfc('%d', sum(trav_wave_all, 2) );
%labels{1} = sprintf('#TW=%d', sum(trav_wave_all(1,:), 2));
xt = get(gca, 'XTick');
text(xt, 8300*ones(size(sum(trav_wave_all, 2)) ), labels, 'FontSize', 20,...
    'HorizontalAlignment', 'center', 'VerticalAlignment','bottom')
text(xt(1)-0.5, 8300, '#TW = ','FontSize', 20,...
    'HorizontalAlignment', 'center', 'VerticalAlignment','bottom');

%}
qsave = 1;
fname = fullfile(save_path_fig, sprintf(...
    't_onset_trav_waves_vs_%s_K12_%d_nruns_%d_digits_%d_boxplot',...
    loopvar_str, K12, nruns, digits));
save_figure(h, 10, 8, fname, '.pdf', qsave);

%% Onset time vs. periodicity
t_data_all = [];
period_data_all = [];

h = figure;
hold on
for i=1:numel(loopvar)
    %i = 1;
    idx = find(period_all(i,:)<Inf);
    scatter(t_onset_all(i, idx), period_all(i, idx) );
    t_data_all = [t_data_all, t_onset_all(i, idx)];
    period_data_all = [period_data_all, period_all(i, idx)];
    
    % calculate correlation
    R = corrcoef(t_onset_all(i, idx), period_all(i, idx));
    fprintf('N = %d, rho = %.2f \n', N_all(i), R(1,2) );
end
legend( sprintfc('N=%d', N_all), 'Location', 'eo' );
% average corr coef
R = corrcoef(t_data_all, period_data_all);
fprintf('Average rho = %.2f \n', R(1,2) )

xlabel('$$t_{onset}$$');
ylabel('Period');
set(gca, 'FontSize', 24);
set(gca, 'XScale', 'log', 'YScale', 'log');
title(sprintf('$\\rho = %.2f$', R(1,2)))
 
qsave = 0;
fname = fullfile(save_path_fig, sprintf(...
    'corr_t_onset_period_var_%s_K12_%d_nruns_%d_digits_%d_',...
    loopvar_str, K12, nruns, digits));
save_figure(h, 10, 8, fname, '.pdf', qsave);