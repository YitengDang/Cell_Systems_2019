%% Analyze trajectory output times and periods
clear all
close all
set(0, 'defaulttextinterpreter', 'latex');
%% Load data
path = 'L:\BN\HY\Shared\Yiteng\two_signals\parameter set 2b';
N = [5 8 10 12 13 15 20].^2;
folders = sprintfc('N%d',N);
%folder = 'N225'; % 'N64', 'N100', 'N225';
%fname_str = sprintf('parameter_set_2b_%s', folder);
fname_str = 'parameter_set_2b';
save_path_fig = 'H:\My Documents\Multicellular automaton\figures\two_signals\analyze_trajectories_t_out_periods';

nN = numel(N);
nruns = 100;
Ilist = 0:0.1:0.7;
I_ini = zeros(nN, numel(Ilist), nruns, 2);
t_out = zeros(nN, numel(Ilist), nruns);
periods = zeros(nN, numel(Ilist), nruns);

%names = {};
for idxN=1:nN
    gz = sqrt(N(idxN));
    dist = init_dist_hex(gz, gz);
    a0 = 1.5;
    folder = folders{idxN};
    listing = dir(fullfile(path, folder));
    num_files = numel(listing)-2; %first two entries are not useful
    %t_out = [];
    %periods = [];

    d = '(\d+)';
    d1 = '0p(\d)0';
    d2 = '(\d+|Inf)';
    %pattern = sprintf('two_signal_mult_%s_chaotic_state_search_t_out_%s_period_%s-v%s',...
    %    folder, d, d2, d); 
    for idx_I=1:numel(Ilist)
        I_ini_temp = Ilist(idx_I)*ones(2,1);
        pattern = strrep(sprintf('two_signal_mult_%s_initiateI1_I_ini_%.2f_%.2f_t_out_%s_period_%s-v%s',...
            folder, I_ini_temp(1), I_ini_temp(2), d, d2, d), '.', 'p');
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
                    cells_ini = cells_hist{1};
                    I_ini(idxN, idx_I, count+1, 1) = moranI(cells_ini(:, 1), a0*dist);
                    I_ini(idxN, idx_I, count+1, 2) = moranI(cells_ini(:, 2), a0*dist);                
                    %t_out(end+1) = str2double(tokens{1}{1});
                    %periods(end+1) = str2double(tokens{1}{2});
                    %I_ini_temp = [str2double(tokens{1}{1}) str2double(tokens{1}{2})]/10;
                    %I_ini(idx_N, idx_I, count(idx_I)+1, :) = I_ini_temp(1);
                    t_out(idxN, idx_I, count+1) = str2double(tokens{1}{1});
                    periods(idxN, idx_I, count+1) = str2double(tokens{1}{2});
                    count = count + 1;
                    %names{count} = name;
                end
            end
        end

    end

end

%% Plot period histogram
% (for single N value at a time)
for idxN=1:nN
    %idxN = 5; % 
    thisN = N(idxN);
    p_data = periods(idxN, :, :);
    uniq = unique(p_data);
    C = categorical(p_data, uniq); 
    h1 = figure(1);
    histogram(C, 'Normalization', 'count');
    xlabel('period');
    ylabel('count');
    title('800 simulations');
    set(gca,'FontSize', 24);
    set(h1, 'Units', 'Inches', 'Position', [1 1 9 8]);

    qsave = 0;
    if qsave
        fname = fullfile(save_path_fig, strcat(fname_str, '_N', num2str(thisN), '_period_hist'));

        save_figure(h1, 9, 8, fname, '.pdf');
    end
    close all
end
%% Plot t_out histogram
% (for single N value at a time)
for idxN=1:nN
    %idxN = 1; % 
    thisN = N(idxN);
    t_data = t_out(idxN, :, :);
    h2 = figure(2);
    tmax = ceil(max(t_data(:))/10)*10;
    edges = 0:tmax/20:tmax;
    histogram(t_data, edges, 'Normalization', 'probability');
    xlabel('$$t_{out}$$');
    ylabel('probability');
    xlim([0 tmax]);
    set(gca,'FontSize', 24);
    set(h2, 'Units', 'Inches', 'Position', [1 1 9 8]);

    qsave = 1;
    if qsave
        fname = fullfile(save_path_fig, strcat(fname_str, '_N', num2str(thisN), '_t_out_hist'));
        save_figure(h2, 9, 8, fname, '.pdf');
    end
    close all
end

%% Scatter t_out against N
% establish correlation and fit linear function
h31 = figure(31);
hold on
N_data = zeros(size(t_data));
for idxN=1:nN
    N_data(idxN, :, :) = N(idxN)*ones(1, size(t_out, 2)*size(t_out, 3));
    t_data = t_out(idxN, :);
end
x = N_data(:);
y = t_out(:);
X = [ones(length(x),1) x];
b = X\y;
yFit = X*b;
Rsqr =  1 - sum((y - yFit).^2)/sum((y - mean(y)).^2);
scatter(x, y);
plot(x, yFit);
[corr_N_t_out,  pval] = corr(N_data(:), t_out(:));
fprintf('Corr(N, t_out) = %.2f, pval = %.2f \n', corr_N_t_out, pval)
xlabel('N');
ylabel('$$t_{out}$$');
ylim([0 14000]);
set(gca,'FontSize', 24);
set(h31, 'Units', 'Inches', 'Position', [1 1 9 8]);
title(sprintf('$$\\rho = %.2f, b = %.2f, R^2 = %.2f$$', corr_N_t_out, b(2), Rsqr), ...
    'FontSize', 24);

qsave = 0;
if qsave
    fname = fullfile(save_path_fig, strcat(fname_str, '_scatter_N_vs_t_out'));
    save_figure(h31, 9, 8, fname, '.pdf');
end
%% correlation for each I
N_data_by_I = zeros(numel(Ilist), size(t_out, 1), size(t_out, 3));
t_out_by_I = zeros(numel(Ilist), size(t_out, 1), size(t_out, 3));
for i=1:numel(Ilist)
    N_data_by_I(i, :) = repmat(N, 1, 100); 
    t_out_by_I(i, :, :) = t_out(:, i, :);
end
[corr_N_t_out_by_I,  pvals] = corr(N_data_by_I(1, :)', t_out_by_I(:, :)');
fprintf('I | Corr(N, t_out) | pval \n');
for i=1:numel(Ilist)
    fprintf('%.2f | %.2f | %.2f \n', Ilist(i), corr_N_t_out_by_I(i), pvals(i))
end
%% Plot <t_out> against N 

%N = [5 8 10 15].^2;
N = [5 8 10 12 13 15 20].^2;

h3 = figure(3);
hold on
%scatter(N, mean(t_out(1:6,:), 2), 100, 'b');
errorbar(N, mean(t_out(:,:), 2), std(t_out(:,:), 1, 2), 'bx');

yFit = [ones(length(N), 1) N']*b;
plot(N, yFit, 'r-');

xlabel('$$N$$');
ylabel('$$\langle t_{out} \rangle$$');
set(gca,'FontSize', 24);
set(h3, 'Units', 'Inches', 'Position', [1 1 9 8]);
title(sprintf('$$\\rho = %.2f, R^2 = %.2f$$', corr_N_t_out, Rsqr), ...
    'FontSize', 24);

qsave = 1;
if qsave
    fname = fullfile(save_path_fig, strcat(fname_str, '_t_out_avg_vs_N_with_fit'));
    save_figure(h3, 9, 8, fname, '.pdf');
end


%% plot periods against N
% Required: load multiple N
N_data = repmat(N', 1, numel(Ilist)*nruns);
periods_stat = zeros(nN, 2);
for i=1:nN
    idx = find(periods(i, :)~=Inf);
    periods_stat(i, 1) = mean(periods(i, idx));
    periods_stat(i, 2) = std(periods(i, idx));
end

x = N_data(:);
y = periods(:);
[corr_N_period,  pval] = corr(x(y~=Inf), y(y~=Inf));
fprintf('Corr(N, period) = %.2f, pval = %.2f \n', corr_N_period, pval)

h4 = figure(4);
hold on

scatter(x(:), y(:), 100);
plot(N, periods_stat(:, 1), 'rd', 'MarkerSize', 20);
%errorbar(N, periods_stat(:, 1), periods_stat(:, 2), 'ro');
set(gca, 'YScale', 'log');
xlabel('$$N$$');
ylabel('Period');
set(gca,'FontSize', 24);
set(h4, 'Units', 'Inches', 'Position', [1 1 9 8]);

qsave = 1;
if qsave
    fname = fullfile(save_path_fig, strcat(fname_str, '_scatter_N_vs_periods'));
    save_figure(h4, 9, 8, fname, '.pdf');
end

%% Plot fraction periodic against N and I
% against N
f_periodic = sum(periods~=Inf, 3)/nruns;
h51 = figure(51);
hold on
scatter(N, mean(f_periodic, 2), 100, 'filled');
xlabel('$$N$$');
ylabel('Fraction periodic');
set(gca,'FontSize', 24);
set(h51, 'Units', 'Inches', 'Position', [1 1 9 8]);
ylim([0 1]);

qsave = 1;
if qsave
    fname = fullfile(save_path_fig, strcat(fname_str, '_frac_periodic_vs_N'));
    save_figure(h51, 9, 8, fname, '.pdf');
end
%% against I
h52 = figure(52);
hold on
scatter(Ilist, mean(f_periodic, 1), 100, 'filled');
xlabel('$$I_{ini}$$');
ylabel('Fraction periodic');
set(gca,'FontSize', 24);
set(h52, 'Units', 'Inches', 'Position', [1 1 9 8]);
%xlim([0 0.701]);
ylim([0 1]);
qsave = 1;
if qsave
    fname = fullfile(save_path_fig, strcat(fname_str, '_frac_periodic_vs_I'));
    save_figure(h52, 9, 8, fname, '.pdf');
end
%% against N and I
h53 = figure(5);
hold on

fdata = sum(periods~=Inf, 3)/nruns;

for idxN=1:nN
    plot(Ilist, fdata(idxN, :), 'o--', 'MarkerSize', 8, 'LineWidth', 1.5);
end
xlabel('$$I_{ini}$$');
ylabel('Fraction periodic');
set(gca,'FontSize', 24);
set(h53, 'Units', 'Inches', 'Position', [1 1 9 8]);
legend(folders, 'Location', 'se');
%xlim([0 0.701]);

qsave = 1;
if qsave
    fname = fullfile(save_path_fig, strcat(fname_str, '_frac_periodic_vs_I_N'));
    save_figure(h53, 9, 8, fname, '.pdf');
end

%% Scatter t_out vs I
h61 = figure(61);
hold on
corr_I_t_out = zeros(nN, 2); % correlations between I1/I2 and t_out
pvals = zeros(nN, 2);
fprintf('Corr(I1, t_out) | Corr(I2, t_out) | pval(I1,t_out) | pval(I2, t_out) \n')
for idxN=1:nN
    I1 = I_ini(idxN, :, :, 1);
    I2 = I_ini(idxN, :, :, 2);
    t_data = t_out(idxN, :);
    scatter3(I1(:), I2(:), t_data(:));
    [corr_I_t_out(idxN, 1), pvals(idxN, 1)] = corr(I1(:), t_data(:)); 
    [corr_I_t_out(idxN, 2), pvals(idxN, 2)] = corr(I2(:), t_data(:));
    fprintf('N=%d | %.2f | %.2f | %.2f | %.2f \n', N(idxN), ...
        corr_I_t_out(idxN, 1), corr_I_t_out(idxN, 2), pvals(idxN, 1), pvals(idxN, 2))
end
xlabel('$$I_1(t=0)$$');
ylabel('$$I_2(t=0)$$');
zlabel('$$t_{out}$$');
set(gca,'FontSize', 24);
set(h61, 'Units', 'Inches', 'Position', [1 1 9 8]);
legend(folders, 'Location', 'east');
xlim([0 1]);
view(45, 30);
text(0.15, 0.15, 18000, ...
    {sprintf('N=%d: $$\\rho$$ = (%.2f, %.2f), p = (%.2f, %.2f)', N(1),...
    corr_I_t_out(1,1), corr_I_t_out(1,2), pvals(1, 1), pvals(1, 2)), ...
    sprintf('N=%d: $$\\rho$$ = (%.2f, %.2f), p = (%.2f, %.2f)', N(2),...
    corr_I_t_out(2,1), corr_I_t_out(2,2), pvals(2, 1), pvals(2, 2)),...
    sprintf('N=%d: $$\\rho$$ = (%.2f, %.2f), p = (%.2f, %.2f)', N(3),...
    corr_I_t_out(3,1), corr_I_t_out(3,2), pvals(3, 1), pvals(3, 2)),...
    sprintf('N=%d: $$\\rho$$ = (%.2f, %.2f), p = (%.2f, %.2f)', N(4),...
    corr_I_t_out(4,1), corr_I_t_out(4,2), pvals(4, 1), pvals(4, 2)),...
    sprintf('N=%d: $$\\rho$$ = (%.2f, %.2f), p = (%.2f, %.2f)', N(5),...
    corr_I_t_out(5,1), corr_I_t_out(5,2), pvals(5, 1), pvals(5, 2)),...
    sprintf('N=%d: $$\\rho$$ = (%.2f, %.2f), p = (%.2f, %.2f)', N(6),...
    corr_I_t_out(6,1), corr_I_t_out(6,2), pvals(6, 1), pvals(6, 2)),...
    sprintf('N=%d: $$\\rho$$ = (%.2f, %.2f), p = (%.2f, %.2f)', N(7),...
    corr_I_t_out(7,1), corr_I_t_out(7,2), pvals(7, 1), pvals(7, 2)),...
    }, 'FontSize', 16)

qsave = 1;
if qsave
    fname = fullfile(save_path_fig, strcat(fname_str, '_scatter_I_ini_vs_t_out'));
    save_figure(h61, 9, 8, fname, '.pdf');
end

%% Plot <t_out> vs I
h6 = figure(6);
hold on

t_mean = mean(t_out, 3);  
t_std = std(t_out, 1, 3);
for idxN=1:nN
    %plot(Ilist, t_mean(idxN, :), 'o--', 'MarkerSize', 8, 'LineWidth', 1.5);
    errorbar(Ilist, t_mean(idxN, :), t_std(idxN, :), 'o--');
end
xlabel('$$I_{ini}$$');
ylabel('$$\langle t_{out} \rangle$$');
set(gca,'FontSize', 24);
set(h6, 'Units', 'Inches', 'Position', [1 1 9 8]);
legend(folders, 'Location', 'eastoutside');
xlim([0 0.8]);

qsave = 1;
if qsave
    fname = fullfile(save_path_fig, strcat(fname_str, '_t_out_avg_vs_I_ini_ebars'));
    save_figure(h6, 9, 8, fname, '.pdf');
end

%% Plot periods vs I
h7 = figure(7);
hold on
corr_I_period = zeros(nN, 2); % correlations between I1/I2 and t_out
pvals = zeros(nN, 2);
fprintf('Corr | (I1, period) | (I2, period) \n')
for idxN=1:nN
    I1 = squeeze(I_ini(idxN, :, :, 1));
    I2 = squeeze(I_ini(idxN, :, :, 2));
    p_data = periods(idxN, :);
    idx = p_data~=Inf;
    scatter3(I1(idx), I2(idx), p_data(idx));
    [corr_I_period(idxN, 1), pvals(idxN, 1)] = corr(I1(idx)', p_data(idx)'); 
    [corr_I_period(idxN, 2), pvals(idxN, 2)] = corr(I2(idx)', p_data(idx)');
    fprintf('N=%d | %.2f | %.2f \n', N(idxN), corr_I_period(idxN, 1), corr_I_period(idxN, 2))
end
xlabel('$$I_1(t=0)$$');
ylabel('$$I_2(t=0)$$');
zlabel('period');
set(gca,'FontSize', 24);
set(h7, 'Units', 'Inches', 'Position', [1 1 9 8]);
legend(folders, 'Location', 'e');
xlim([0 1]);
view(45, 30);
text(0.15, 0.15, 6000, ...
    {sprintf('N=%d: $$\\rho$$ = (%.2f, %.2f), p = (%.2f, %.2f)', N(1),...
    corr_I_period(1,1), corr_I_period(1,2), pvals(1, 1), pvals(1, 2)), ...
    sprintf('N=%d: $$\\rho$$ = (%.2f, %.2f), p = (%.2f, %.2f)', N(2),...
    corr_I_period(2,1), corr_I_period(2,2), pvals(2, 1), pvals(2, 2)),...
    sprintf('N=%d: $$\\rho$$ = (%.2f, %.2f), p = (%.2f, %.2f)', N(3),...
    corr_I_period(3,1), corr_I_period(3,2), pvals(3, 1), pvals(3, 2)),...
    sprintf('N=%d: $$\\rho$$ = (%.2f, %.2f), p = (%.2f, %.2f)', N(4),...
    corr_I_period(4,1), corr_I_period(4,2), pvals(4, 1), pvals(4, 2)),...
    sprintf('N=%d: $$\\rho$$ = (%.2f, %.2f), p = (%.2f, %.2f)', N(5),...
    corr_I_period(5,1), corr_I_period(5,2), pvals(5, 1), pvals(5, 2)),...
    sprintf('N=%d: $$\\rho$$ = (%.2f, %.2f), p = (%.2f, %.2f)', N(6),...
    corr_I_period(6,1), corr_I_period(6,2), pvals(6, 1), pvals(6, 2)),...
    sprintf('N=%d: $$\\rho$$ = (%.2f, %.2f), p = (%.2f, %.2f)', N(7),...
    corr_I_period(7,1), corr_I_period(7,2), pvals(7, 1), pvals(7, 2)),...
    }, 'FontSize', 16)

qsave = 1;
if qsave
    fname = fullfile(save_path_fig, strcat(fname_str, '_scatter_I_ini_vs_period'));
    save_figure(h7, 9, 8, fname, '.pdf');
end
