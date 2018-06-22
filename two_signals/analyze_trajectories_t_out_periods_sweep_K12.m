%% Analyze trajectory output times and periods
clear all
close all
set(0, 'defaulttextinterpreter', 'latex');
%% Load files
% load path
%path = 'L:\BN\HY\Shared\Yiteng\two_signals\parameter set 2b';
path = 'L:\BN\HY\Shared\Yiteng\two_signals\sweep K12';

% save path & filename
fname_str = 'parameter_set_2b';
save_path_fig = fullfile('H:\My Documents\Multicellular automaton\figures\two_signals',...
    'analyze_trajectories_vs_K12');

% Parameters
%N = [5 8 10 12 13 15 20].^2;
N = 225;
nN = numel(N);

K12_all = 5:1:13;
nK12 = numel(K12_all);
a0 = 1.5;
nruns = 500;
Ilist = 0:0.1:0.7;
%I_ini = zeros(nN, numel(Ilist), nruns, 2);

folders = sprintfc('N%d',N);
%folder = 'N225'; % 'N64', 'N100', 'N225';
%fname_str = sprintf('parameter_set_2b_%s', folder);

% Load and store data
t_out = zeros(nK12, nruns);
periods = zeros(nK12, nruns);

%names = {};
for idxK12=1:nK12
for idxN=1:nN
    K12 = K12_all(idxK12);
    gz = sqrt(N(idxN));
    dist = init_dist_hex(gz, gz);
    
    %folder = folders{idxN};
    folder = '';
    listing = dir(fullfile(path, folder));
    num_files = numel(listing)-2; %first two entries are not useful
    %t_out = [];
    %periods = [];

    d = '(\d+)';
    d1 = '0p(\d)0';
    d2 = '(\d+|Inf)';
    %pattern = sprintf('two_signal_mult_%s_chaotic_state_search_t_out_%s_period_%s-v%s',...
    %    folder, d, d2, d); 
    %for idx_I=1:numel(Ilist)
        %I_ini_temp = Ilist(idx_I)*ones(2,1);
        %pattern = strrep(sprintf('two_signal_mult_%s_initiateI1_I_ini_%.2f_%.2f_t_out_%s_period_%s-v%s',...
        %    folder, I_ini_temp(1), I_ini_temp(2), d, d2, d), '.', 'p');
        pattern = strrep(sprintf('two_signal_mult_N%d_initiateI0_K12_%d_t_out_%s_period_%s-v%s',...
            N, K12, d, d2, d), '.', 'p');
        disp(pattern);
        count = 0;
        for i = 1:num_files
            filename = listing(i+2).name;
            % remove extension and do not include txt files
            [~,name,ext] = fileparts(filename);
            if strcmp(ext, '.mat')
                [~, tokens] = regexp(name, pattern, 'match', 'tokens');
                if ~isempty(tokens)
                    load(fullfile(path, folder, name), 'cells_hist');
                    disp(name);
                    cells_ini = cells_hist{1};
                    %I_ini(idxN, idx_I, count+1, 1) = moranI(cells_ini(:, 1), a0*dist);
                    %I_ini(idxN, idx_I, count+1, 2) = moranI(cells_ini(:, 2), a0*dist);                
                    %t_out(end+1) = str2double(tokens{1}{1});
                    %periods(end+1) = str2double(tokens{1}{2});
                    %I_ini_temp = [str2double(tokens{1}{1}) str2double(tokens{1}{2})]/10;
                    %I_ini(idx_N, idx_I, count(idx_I)+1, :) = I_ini_temp(1);
                    t_out(idxK12, count+1) = str2double(tokens{1}{1});
                    periods(idxK12, count+1) = str2double(tokens{1}{2});
                    count = count + 1;
                    %names{count} = name;
                end
            end
        end

    %end

end
end

%% Plot t_out histogram
% (for single N value at a time)
%for idxN=1:nN
    %idxN = 1; % 
    %thisN = N(idxN);
    %t_data = t_out(idxN, :, :);

% Plot together
h2 = figure(2);
hold on
%tmax = ceil(max(t_data(:))/10)*10;
%edges = 0:tmax/20:tmax;
nbins = 10;
for i=1:numel(K12_all)
    histogram(t_out(i,:), nbins, 'Normalization', 'probability');
end
xlabel('$$t_{out}$$');
ylabel('probability');
%xlim([0 tmax]);
set(gca,'FontSize', 24);
set(gca, 'XScale', 'log');
set(h2, 'Units', 'Inches', 'Position', [1 1 9 8]);
legend(sprintfc('$$K_{12}=%d$$',K12_all), 'Interpreter', 'latex');

qsave = 0;
if qsave
    fname = fullfile(save_path_fig, strcat(fname_str, '_t_out_hist_all_',...
        num2str(nruns) ,'runs'));
    save_figure(h2, 9, 8, fname, '.pdf');
end
    %close all
%end

%% Scatter t_out against K12
% establish correlation and fit linear function
h31 = figure(3);
hold on
K12_data = zeros(size(t_out));
for idxK12=1:nK12
    K12_data(idxK12, :, :) = K12_all(idxK12)*ones(1, size(t_out, 2));
end
x = K12_data(:);
y = t_out(:);

% fit linear function
%X = [ones(length(x),1) x];
%b = X\y;
%yFit = X*b;
%Rsqr =  1 - sum((y - yFit).^2)/sum((y - mean(y)).^2);
%plot(x, yFit);

% plot data
scatter(x, y);
%[corr_N_t_out,  pval] = corr(K12_data(:), t_out(:));
%fprintf('Corr(N, t_out) = %.2f, pval = %.2f \n', corr_N_t_out, pval)
xlabel('$$K_{12}$$');
ylabel('$$t_{out}$$');
%ylim([0 14000]);
set(gca,'FontSize', 24);
set(gca, 'YScale', 'log');
set(h31, 'Units', 'Inches', 'Position', [1 1 9 8]);
%title(sprintf('$$\\rho = %.2f, b = %.2f, R^2 = %.2f$$', corr_N_t_out, b(2), Rsqr), ...
%    'FontSize', 24);

qsave = 0;
if qsave
    fname = fullfile(save_path_fig, strcat(fname_str, '_K12_vs_t_out_scatter_', ...
        num2str(nruns) ,'runs'));
    save_figure(h31, 9, 8, fname, '.pdf');
end

% Plot <t_out> against K12
%N = [5 8 10 15].^2;
%N = [5 8 10 12 13 15 20].^2;

h3 = figure(3);
hold on
%scatter(N, mean(t_out(1:6,:), 2), 100, 'b');
%errorbar(K12_all, mean(t_out, 2), std(t_out, 1, 2), 'rd', 'LineWidth', 2, 'MarkerSize', 20);
plot(K12_all, mean(t_out, 2), 'rd', 'LineWidth', 2, 'MarkerSize', 20);

%yFit = [ones(length(N), 1) N']*b;
%plot(N, yFit, 'r-');

xlabel('$$K_{12}$$');
%ylabel('$$\langle t_{out} \rangle$$');
set(gca,'FontSize', 24);
%set(gca, 'YScale', 'log');
set(h3, 'Units', 'Inches', 'Position', [1 1 9 8]);
%title(sprintf('$$\\rho = %.2f, R^2 = %.2f$$', corr_N_t_out, Rsqr), ...
%    'FontSize', 24);
ylim(10.^[0 5]);
set(gca, 'YTick', 10.^[0:5]);

qsave = 1;
if qsave
    fname = fullfile(save_path_fig, strcat(fname_str, '_K12_vs_t_out_scatter_errorbar_', ...
        num2str(nruns) ,'runs'));
    save_figure(h3, 9, 8, fname, '.pdf');
end

%% Plot <t_out> against K12 - boxplot

h32 = figure(32);
hold on
%scatter(repmat(1:numel(K12_all), 1, nruns), y);
boxplot(t_out', K12_all); %, 'PlotStyle', 'compact' );

xlabel('$$K_{12}$$');
ylabel('$$t_{out}$$');
set(gca, 'YScale', 'log');

set(gca,'FontSize', 24);
set(h32, 'Units', 'Inches', 'Position', [1 1 9 8]);

ylim(10.^[0 5]);
set(gca, 'YTick', 10.^[0:5]);

qsave = 1;
if qsave
    fname = fullfile(save_path_fig, strcat(fname_str, '_t_out_vs_K12_box_plot_', ...
        num2str(nruns) ,'runs', '_log'));
    save_figure(h32, 9, 8, fname, '.pdf');
end

%% Plot fraction of periodic trajectories vs K12
% and fraction with period 4
f_nonperiodic = sum(periods==Inf, 2)/nruns;
f_period4 = sum(periods==4, 2)/nruns;

h5 = figure(5);
hold on
yyaxis left
plot(K12_all, f_nonperiodic, 'b--', 'MarkerSize', 10);
scatter(K12_all, f_nonperiodic, 100, 'b', 'filled');

yyaxis right
plot(K12_all, f_period4, 'r--', 'MarkerSize', 10);
scatter(K12_all, f_period4, 100, 'r', 'filled');
%scatter(K12_all, f_nonperiodic, 100, 'filled');

yyaxis left
xlabel('$$K_{12}$$');
ylabel('Fraction non-periodic');
ylim([0 1.05]);

yyaxis right
ylabel('Fraction period 4');
ylim([0 1.05]);

set(gca,'FontSize', 24);
set(h5, 'Units', 'Inches', 'Position', [1 1 9 8]);

qsave = 1;
if qsave
    fname = fullfile(save_path_fig, strcat(fname_str, '_frac_periodic_vs_K12_',...
        num2str(nruns), 'runs_v2b'));
    save_figure(h5, 9, 8, fname, '.pdf');
end
%% Plot period histogram
% (for single N value at a time)
for idxK12=1:nK12
    %idxN = 5; % 
    %thisN = N(idxN);
    
    p_data = periods(idxK12, :);
    uniq = unique(p_data);
    C = categorical(p_data, uniq); 
    h1 = figure(idxK12);
    histogram(C, 'Normalization', 'count');
    xlabel('period');
    ylabel('count');
    title(sprintf('%d simulations', nruns));
    set(gca,'FontSize', 24);
    set(h1, 'Units', 'Inches', 'Position', [1 1 9 8]);

    qsave = 0;
    if qsave
        fname = fullfile(save_path_fig, strcat(fname_str, '_K12',...
            num2str(K12_all(idxK12)), '_period_hist'));
        save_figure(h1, 9, 8, fname, '.pdf');
    end
    %close all
end

%% plot periods against K12
K12_data = repmat(K12_all', 1, nruns);
periods_stat = zeros(nK12, 2);
for i=1:nK12
    idx = find(periods(i, :)~=Inf);
    periods_stat(i, 1) = mean(periods(i, idx));
    periods_stat(i, 2) = std(periods(i, idx));
end

x = K12_data(:);
y = periods(:);
%[corr_N_period,  pval] = corr(x(y~=Inf), y(y~=Inf));
%fprintf('Corr(N, period) = %.2f, pval = %.2f \n', corr_N_period, pval)

h4 = figure(4);
hold on
scatter(x(:), y(:), 100);
plot(K12_all, periods_stat(:, 1), 'rd', 'MarkerSize', 20);
%errorbar(N, periods_stat(:, 1), periods_stat(:, 2), 'ro');
%set(gca, 'YScale', 'log');
xlim([K12_all(1)-1 K12_all(end)+1]);
xlabel('$$K_{12}$$');
ylabel('Period');
set(gca,'FontSize', 24);
set(gca,'YScale','log');
set(h4, 'Units', 'Inches', 'Position', [1 1 9 8]);

qsave = 1;
if qsave
    fname = fullfile(save_path_fig, strcat(fname_str, '_K12_vs_periods_scatter_mean_',...
        num2str(nruns), 'runs'));
    save_figure(h4, 9, 8, fname, '.pdf');
end