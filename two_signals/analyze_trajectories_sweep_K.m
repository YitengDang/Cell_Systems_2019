%% Analyze trajectory output times and periods
clear all
close all
set(0, 'defaulttextinterpreter', 'latex');
%% Load files
% load path
%path = 'L:\BN\HY\Shared\Yiteng\two_signals\parameter set 2b';
%path = 'L:\BN\HY\Shared\Yiteng\two_signals\sweep K12';
path = 'L:\BN\HY\Shared\Yiteng\two_signals\sweep K22 new lattice';
%path = fullfile('L:\BN\HY\Shared\Yiteng\two_signals\sweep K12 new lattice', sprintf('N%d', N)));
%path = fullfile('L:\BN\HY\Shared\Yiteng\two_signals\sweep K12 new lattice', 'N225 mcsteps 1000');
mcsteps = 0;
tmax = 10^5;
lattice_lbl = sprintf('randpos_mcsteps%d', mcsteps);

% save path
%save_path_fig = fullfile('H:\My Documents\Multicellular automaton\figures\two_signals',...
%    'analyze_trajectories_vs_K12');
save_path_fig = fullfile('H:\My Documents\Multicellular automaton\latex\10_two_signals\analyze_trajectories', ...
    'analyze_trajectories_temp');

% Parameters
%N = [5 8 10 12 13 15 20].^2;
N = 625; %225;
nN = numel(N);

K_all = [5:10 15:5:25]; %[7:11 15 20:24];
nK = numel(K_all);
a0 = 0.5;
nruns = 100;
%Ilist = 0:0.1:0.7;
%I_ini = zeros(nN, numel(Ilist), nruns, 2);

folders = sprintfc('N%d',N);
%folder = 'N225'; % 'N64', 'N100', 'N225';
%fname_str = sprintf('parameter_set_2b_%s', folder);

% Load and store data
t_out = zeros(nK, nruns);
periods = zeros(nK, nruns);
n_trav_wave = zeros(nK, 1);
p_av_all = zeros(nK, nruns, 2);
pij_av_all = zeros(nK, nruns, 4);
I_av_all = zeros(nK, nruns, 2);

%names = {};
for idxK=1:nK
    for idxN=1:nN
        K22 = K_all(idxK);
        gz = sqrt(N(idxN));
        dist = init_dist_hex(gz, gz);
        
        fprintf('gz = %d, K_22 = %.2f \n', gz, K22);
        
        %path = fullfile('L:\BN\HY\Shared\Yiteng\two_signals\sweep K12 new lattice', sprintf('N%d', gz^2));

        %folder = folders{idxN};
        %folder = '';
        folder = 'N625 strong int a0 0p5';
        listing = dir(fullfile(path, folder));
        num_files = numel(listing)-2; %first two entries are not useful
        %t_out = [];
        %periods = [];

        t1 = '(\d+)';
        t2 = '0p(\d)0';
        t3 = '(\d+|Inf)';
        t4 = '\w*';
        

        %pattern = sprintf('two_signal_mult_%s_chaotic_state_search_t_out_%s_period_%s-v%s',...
        %    folder, d, d2, d); 
        %for idx_I=1:numel(Ilist)
            %I_ini_temp = Ilist(idx_I)*ones(2,1);
            %pattern = strrep(sprintf('two_signal_mult_%s_initiateI1_I_ini_%.2f_%.2f_t_out_%s_period_%s-v%s',...
            %    folder, I_ini_temp(1), I_ini_temp(2), d, d2, d), '.', 'p');
            %pattern = strrep(sprintf(...
            %     'two_signal_mult_N%d_initiateI0_K12_%d_t_out_%s_period_%s-v%s',...
            %     N, K22, t1, t3, t1), '.', 'p');
                       
            %pattern = strrep(sprintf(...
            %     'two_signal_mult_N%d_initiateI0_randpos_mcsteps0_K12_%d_t_out_%s_period_%s-v%s',...
            %     N, K12, t1, t3, t1), '.', 'p');
            %pattern = strrep(sprintf(...
            %    'two_signal_mult_N%d_initiateI0_randpos_mcsteps%d_K12_%d_t_out_%s_period_%s%s-v%s',...
            %    N, mcsteps, K12, t1, t3, t4, t1), '.', 'p');
            pattern = strrep(sprintf(...
                'two_signal_mult_N%d_initiateI0_randpos_mcsteps%d_K22_%d_t_out_%s_period_%s%s_temp-v%s',...
                N, mcsteps, K22, t1, t3, t4, t1), '.', 'p');
            
            %disp(pattern);
            count = 0;
            for i = 1:num_files
                filename = listing(i+2).name;
                % remove extension and do not include txt files
                [~,name,ext] = fileparts(filename);
                if strcmp(ext, '.mat')
                    [~, tokens] = regexp(name, pattern, 'match', 'tokens');
                    if ~isempty(tokens)
                        % Load data
                        load(fullfile(path, folder, name), 'cells_hist');
                        disp(name);
                        this_t_out = str2double(tokens{1}{1});
                        this_period = str2double(tokens{1}{2});
                        
                        % Store data
                        %cells_ini = cells_hist{1};
                        %I_ini(idxN, idx_I, count+1, 1) = moranI(cells_ini(:, 1), a0*dist);
                        %I_ini(idxN, idx_I, count+1, 2) = moranI(cells_ini(:, 2), a0*dist);                
                        %t_out(end+1) = str2double(tokens{1}{1});
                        %periods(end+1) = str2double(tokens{1}{2});
                        %I_ini_temp = [str2double(tokens{1}{1}) str2double(tokens{1}{2})]/10;
                        %I_ini(idx_N, idx_I, count(idx_I)+1, :) = I_ini_temp(1);
                        t_out(idxK, count+1) = this_t_out;
                        periods(idxK, count+1) = this_period;
                        count = count + 1;
                        %names{count} = name;
                        
                        trav_wave = travelling_wave_test(cells_hist, a0, this_period, this_t_out);
                        n_trav_wave(idxK) = n_trav_wave(idxK) + trav_wave;
                        
                        [p_av, pij_av, I_av] = final_pI_av(cells_hist, a0, this_period, this_t_out);
                        p_av_all(idxK, count+1, :) = p_av;
                        pij_av_all(idxK, count+1, :) = pij_av;
                        I_av_all(idxK, count+1, :) = I_av;
                    end
                end
            end

        %end

    end
end

%% Save data
%folder = 'L:\BN\HY\Shared\Yiteng\two_signals\sweep K12\analyzed_data';
%folder = fullfile('L:\BN\HY\Shared\Yiteng\two_signals\sweep K12 new lattice', 'analyzed_data');
folder = 'L:\BN\HY\Shared\Yiteng\two_signals\sweep K22 new lattice\analyzed_data';
fname_str = sprintf('Sweep_K22_%d_to_%d_N%d_%s', K_all(1), K_all(end), N, lattice_lbl);
fname = fullfile(folder, strcat(fname_str, '_data'));
%save(fname, 'N', 'a0', 'K_all', 'nruns', 't_out', 'periods',...
%    'n_trav_wave', 'p_av_all', 'pij_av_all', 'I_av_all');


%% Plot trajectory types vs K12
% (1) fraction of stationary/non-periodic trajectories
% (2) fraction of travelling waves
% (3) fraction of trajectories with period 4
% (4) fraction of other periodic trajectories

f_stationary = sum(periods==Inf, 2)/nruns;
f_trav_wave = n_trav_wave/nruns;
f_period4 = sum(periods==4, 2)/nruns;
f_period_other = 1 - f_stationary - f_trav_wave - f_period4; 

h = figure(1);
hold on
plot(K_all, f_stationary, 'k--', 'MarkerSize', 10);
scatter(K_all, f_stationary, 100, 'k', 'filled', 'MarkerFaceAlpha', 0.5);

plot(K_all, f_trav_wave, 'b--', 'MarkerSize', 10);
scatter(K_all, f_trav_wave, 100, 'b', 'filled', 'MarkerFaceAlpha', 0.5);

plot(K_all, f_period4, 'r--', 'MarkerSize', 10);
scatter(K_all, f_period4, 100, 'r', 'filled', 'MarkerFaceAlpha', 0.5);

plot(K_all, f_period_other, 'g--', 'MarkerSize', 10);
scatter(K_all, f_period_other, 100, 'g', 'filled', 'MarkerFaceAlpha', 0.5);

% settings
xlim([K_all(1) K_all(end)])
%xlabel('$$K_{12}$$');
xlabel('$$K_{22}$$');
ylabel('Fraction');
set(gca,'FontSize', 24);
set(h, 'Units', 'Inches', 'Position', [1 1 9 8]);
set(gca, 'XTick', K_all(1):4:K_all(end));

qsave = 1;
fname = fullfile(save_path_fig, strcat(fname_str, '_frac_all_types_vs_K_v1_',...
    num2str(nruns), 'runs'));
save_figure(h, 9, 8, fname, '.pdf', qsave);

%% Plot trajectory types vs K12, v2
% (1) trajectories that reach tmax
% (2) trajectories that end up in a stationary state at t<tmax
% (3) trajectories that end up in oscillatory states

f_tmax = sum(t_out == tmax, 2)/nruns;
f_stationary = sum(periods==Inf & t_out<tmax, 2)/nruns;
f_period4 = sum(periods==4, 2)/nruns;
f_period_other = 1 - f_stationary - f_tmax - f_period4; 
%f_osc = sum(periods<Inf, 2)/nruns;

h = figure(1);
hold on
plot(K_all, f_tmax, 'm--', 'MarkerSize', 10);
scatter(K_all, f_tmax, 100, 'm', 'filled', 'MarkerFaceAlpha', 0.5);

plot(K_all, f_stationary, 'k--', 'MarkerSize', 10);
scatter(K_all, f_stationary, 100, 'k', 'filled', 'MarkerFaceAlpha', 0.5);

plot(K_all, f_period4, 'r--', 'MarkerSize', 10);
scatter(K_all, f_period4, 100, 'r', 'filled', 'MarkerFaceAlpha', 0.5);

plot(K_all, f_period_other, 'g--', 'MarkerSize', 10);
scatter(K_all, f_period_other, 100, 'g', 'filled', 'MarkerFaceAlpha', 0.5);

% settings
xlim([K_all(1) K_all(end)])
%xlabel('$$K_{12}$$');
xlabel('$$K_{22}$$');
ylabel('Fraction');
set(gca,'FontSize', 24);
set(h, 'Units', 'Inches', 'Position', [1 1 9 8]);
set(gca, 'XTick', K_all(1):4:K_all(end));

qsave = 1;
fname = fullfile(save_path_fig, strcat(fname_str, '_frac_all_types_vs_K_v2_',...
    num2str(nruns), 'runs_v2'));
save_figure(h, 9, 8, fname, '.pdf', qsave);

%% Plot two at the same time
h = figure(11);
hold on

% plot non-periodic
%
yyaxis left
plot(K_all, f_stationary, 'b--', 'MarkerSize', 10);
scatter(K_all, f_stationary, 100, 'b', 'filled');
ylabel('Fraction non-periodic');
ylim([0 1.05]);
%}
% plot period 4
yyaxis right
plot(K_all, f_period4, 'r--', 'MarkerSize', 10);
scatter(K_all, f_period4, 100, 'r', 'filled');
%scatter(K12_all, f_nonperiodic, 100, 'filled');
yyaxis right
ylabel('Fraction period 4');
ylim([0 1.05]);

% plot travelling waves
%{
yyaxis left
plot(K12_all, n_trav_wave/nruns, 'k--', 'MarkerSize', 10);
scatter(K12_all, n_trav_wave/nruns, 100, 'k', 'filled');
%scatter(K12_all, f_nonperiodic, 100, 'filled');
ylabel('Fraction travelling waves');
ylim([0 1.05]);
%}
% general settings
xlim([K_all(1) K_all(end)])
%xlabel('$$K_{12}$$');
xlabel('$$K_{22}$$');
set(gca,'FontSize', 24);
set(h, 'Units', 'Inches', 'Position', [1 1 9 8]);
set(gca, 'XTick', K_all(1):4:K_all(end));

qsave = 0;
if qsave
    fname = fullfile(save_path_fig, strcat(fname_str, '_frac_nonperiodic_period4_vs_K_',...
        num2str(nruns), 'runs'));
    save_figure(h, 9, 8, fname, '.pdf');
end
%% Plot frac. travelling waves and I(final)
f_stationary = sum(periods==Inf, 2)/nruns;
f_period4 = sum(periods==4, 2)/nruns;
I_av = squeeze(mean(I_av_all, 2));

h = figure(2);

% set axes colors
left_color = [0.8 0 0.8];
right_color = [0.1 0.1 .1];
set(h,'defaultAxesColorOrder', [left_color; right_color]);

hold on

% plot travelling waves
%
yyaxis right
plot(K_all, n_trav_wave/nruns, 'k--', 'MarkerSize', 10);
scatter(K_all, n_trav_wave/nruns, 100, 'k', 'filled', 'MarkerFaceAlpha', 0.5);
ylabel('Fraction travelling waves');
ylim([0 1.05]);
%}
% Plot with average I_final
yyaxis left
plot(K_all, I_av(:,1), 'b--'); %, 'MarkerSize', 10, 'MarkerFaceColor', 'b');
plot(K_all, I_av(:,2), 'r--'); %, 'MarkerSize', 10, 'MarkerFaceColor', 'r');
p1=scatter(K_all, I_av(:,1),  100, 'b', 'filled', 'MarkerFaceAlpha', 0.5);
p2=scatter(K_all, I_av(:,2),  100, 'r', 'filled', 'MarkerFaceAlpha', 0.5);
legend([p1 p2], {'1', '2'}, 'Location', 'nw');
ylim([0 1.05]);
ylabel('Final $$\langle I^{(i)}\rangle$$');

% general settings
%xlabel('$$K_{12}$$');
xlabel('$$K_{22}$$');
xlim([K_all(1) K_all(end)])
set(gca,'FontSize', 24);
set(h, 'Units', 'Inches', 'Position', [1 1 9 8]);
set(gca, 'XTick', K_all(1):4:K_all(end));

qsave = 1;
if qsave
    fname = fullfile(save_path_fig, strcat(fname_str, '_frac_trav_wave_and_Ifinal_vs_K_',...
        num2str(nruns), 'runs'));
    save_figure(h, 9, 8, fname, '.pdf');
end

%% Plot average final p
p_av = squeeze(mean(p_av_all, 2));

h = figure(31);
hold on
plot(K_all, p_av(:,1), 'b--'); %, 'MarkerSize', 10, 'MarkerFaceColor', 'b');
plot(K_all, p_av(:,2), 'r--'); %, 'MarkerSize', 10, 'MarkerFaceColor', 'r');
p1=scatter(K_all, p_av(:,1), 100, 'b', 'filled', 'MarkerFaceAlpha', 0.5);
p2=scatter(K_all, p_av(:,2), 100, 'r', 'filled', 'MarkerFaceAlpha', 0.5);
legend([p1 p2], {'1', '2'});

% settigns
xlim([K_all(1) K_all(end)])
ylim([0 1]);
xlabel('$$K_{12}$$');
ylabel('Final $$\langle p^{(i)}\rangle$$');
set(gca,'FontSize', 24);
set(h, 'Units', 'Inches', 'Position', [1 1 9 8]);
set(gca, 'XTick', K_all(1):4:K_all(end));

qsave = 1;
if qsave
    fname = fullfile(save_path_fig, strcat(fname_str, '_p_av_vs_K12_',...
        num2str(nruns), 'runs'));
    save_figure(h, 9, 8, fname, '.pdf');
end

%% Plot average final p_ij
pij_av = mean(pij_av_all, 2);

h = figure(32);
hold on
plot(K_all, pij_av(:,1), 'w--')%, 'MarkerSize', 10, 'MarkerFaceColor', 'b');
plot(K_all, pij_av(:,2), 'y--')%, 'MarkerSize', 10, 'MarkerFaceColor', 'r');
plot(K_all, pij_av(:,3), 'b--')%, 'MarkerSize', 10, 'MarkerFaceColor', 'r');
plot(K_all, pij_av(:,4), 'k--')%, 'MarkerSize', 10, 'MarkerFaceColor', 'r');
p1=scatter(K_all, pij_av(:,1), 100, 'w', 'filled', 'MarkerFaceAlpha', 1, 'MarkerEdgeColor', 'k');
p2=scatter(K_all, pij_av(:,2), 100, 'y', 'filled', 'MarkerFaceAlpha', 0.5, 'MarkerEdgeColor', 'k');
p3=scatter(K_all, pij_av(:,3), 100, 'b', 'filled', 'MarkerFaceAlpha', 0.5, 'MarkerEdgeColor', 'k');
p4=scatter(K_all, pij_av(:,4), 100, 'k', 'filled', 'MarkerFaceAlpha', 0.5, 'MarkerEdgeColor', 'k');

% settings
set(gca, 'Color', [0.8 0.8 0.8]);
legend([p1 p2 p3 p4], {'(0,0)', '(1,0)', '(0,1)', '(1,1)'}, 'Location', 'eo');
xlim([K_all(1) K_all(end)])
xlabel('$$K_{12}$$');
ylabel('Final $$\langle p^{(i,j)}\rangle$$');
set(gca,'FontSize', 24);
set(h, 'Units', 'Inches', 'Position', [1 1 12 8]);
set(gca, 'XTick', K_all(1):4:K_all(end));

qsave = 1;
if qsave
    h.Color = 'white';
    h.InvertHardcopy = 'off';
    fname = fullfile(save_path_fig, strcat(fname_str, '_pij_av_vs_K_',...
        num2str(nruns), 'runs'));
    save_figure(h, 10, 8, fname, '.pdf');
end

%% Scatter t_out against K12
% data processing
K_data = zeros(size(t_out));
for idxK=1:nK
    K_data(idxK, :, :) = K_all(idxK)*ones(1, size(t_out, 2));
end
t_idx = t_out<tmax;
x1 = K_data(t_idx);
y1 = t_out(t_idx);
x2 = K_data(~t_idx);
y2 = t_out(~t_idx);

% establish correlation and fit linear function
% fit linear function
%X = [ones(length(x),1) x];
%b = X\y;
%yFit = X*b;
%Rsqr =  1 - sum((y - yFit).^2)/sum((y - mean(y)).^2);
%plot(x, yFit);

% plot data
h = figure(41);
hold on
scatter(x1, y1);
scatter(x2, y2, 'r');
plot([K_all(1) K_all(end)], [tmax tmax], 'r--');
%[corr_N_t_out,  pval] = corr(K12_data(:), t_out(:));
%fprintf('Corr(N, t_out) = %.2f, pval = %.2f \n', corr_N_t_out, pval)

qsave = 0;
if qsave
    fname = fullfile(save_path_fig, strcat(fname_str, '_K12_vs_t_out_scatter_', ...
        num2str(nruns) ,'runs'));
    save_figure(h, 9, 8, fname, '.pdf');
end

% Plot <t_out> against K12
%N = [5 8 10 15].^2;
%N = [5 8 10 12 13 15 20].^2;

h = figure(41);
hold on
%scatter(N, mean(t_out(1:6,:), 2), 100, 'b');
%errorbar(K12_all, mean(t_out, 2), std(t_out, 1, 2), 'rd', 'LineWidth', 2, 'MarkerSize', 20);
plot(K_all, mean(t_out, 2), 'rd', 'LineWidth', 2, 'MarkerSize', 20);

%yFit = [ones(length(N), 1) N']*b;
%plot(N, yFit, 'r-');

% figure settings
%xlabel('$$K_{12}$$');
xlabel('$$K_{22}$$');
ylabel('$$t_{out}$$');
xlim([K_all(1) K_all(end)]);
ylim(10.^[0 5]);
set(gca, 'YTick', 10.^(0:5));
set(gca,'FontSize', 24);
set(gca, 'YScale', 'log');
set(h, 'Units', 'Inches', 'Position', [1 1 9 8]);
%title(sprintf('$$\\rho = %.2f, b = %.2f, R^2 = %.2f$$', corr_N_t_out, b(2), Rsqr), ...
%    'FontSize', 24);

qsave = 1;
if qsave
    fname = fullfile(save_path_fig, strcat(fname_str, '_K_vs_t_out_mean_scatter_', ...
        num2str(nruns) ,'runs'));
    save_figure(h, 9, 8, fname, '.pdf');
end

%% Plot <t_out> against K12 - boxplot
h = figure(42);
hold on
%scatter(repmat(1:numel(K12_all), 1, nruns), y);
p=boxplot(t_out', K_all); %, 'PlotStyle', 'compact' );

xlabel('$$K_{12}$$');
ylabel('$$t_{out}$$');
set(gca, 'YScale', 'log');

set(gca,'FontSize', 24);
set(h, 'Units', 'Inches', 'Position', [1 1 9 8]);

ylim(10.^[0 5]);
set(gca, 'XTick', 1:2:numel(K_all), 'xticklabel', K_all(1):2:K_all(end));
set(gca, 'YTick', 10.^(0:5));

qsave = 0;
if qsave
    fname = fullfile(save_path_fig, strcat(fname_str, '_t_out_vs_K_box_plot_', ...
        num2str(nruns) ,'runs', '_log'));
    save_figure(h, 9, 8, fname, '.pdf');
end

%% plot periods against K12
K_data = repmat(K_all', 1, nruns);

x = K_data(:); %+rand(size(K12_data(:)))*0.1; % randomize x positions
y = periods(:);
%[corr_N_period,  pval] = corr(x(y~=Inf), y(y~=Inf));
%fprintf('Corr(N, period) = %.2f, pval = %.2f \n', corr_N_period, pval)

periods_stat = zeros(nK, 2);
h = figure(51);
hold on

uniq_periods = unique(periods);
uniq_periods = sort(uniq_periods(1:end-1), 'descend');
uniq_periods_str = string(uniq_periods);
plot_list = zeros(numel(uniq_periods)-1, 1);
for i=1:numel(uniq_periods)
    plot_list(i) = plot([K_all(1)-1 K_all(end)+1], [uniq_periods(i) uniq_periods(i)], '--', 'LineWidth', 0.5);
end
for i=1:nK
    K22 = K_all(i);
    % get mean & std
    idx = find(periods(i, :)~=Inf);
    periods_stat(i, 1) = mean(periods(i, idx));
    periods_stat(i, 2) = std(periods(i, idx));
    
    % count unique periods
    uniq_periods = unique(periods(i, :));
    C = categorical(periods(i, idx));
    [counts, cats] = histcounts(C);
    num_traj = sum(counts);
    
    if ~isempty(cats)
        %fprintf('K12 = %d, \n', K12);
        x = repmat(K22, 1, numel(cats));
        scatter(x, str2double(cats), 6*counts, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'b');
    end
end

legend(plot_list, uniq_periods_str, 'Location', 'eo')
set(gca, 'YScale', 'log');
xlim([K_all(1)-1 K_all(end)+1]);
%xlabel('$$K_{12}$$');
xlabel('$$K_{22}$$');
ylabel('Period');
set(gca,'FontSize', 24);
%set(gca,'YScale','log');
set(h, 'Units', 'Inches', 'Position', [1 1 9 8]);

qsave = 1;
if qsave
    fname = fullfile(save_path_fig, strcat(fname_str, '_K_vs_periods_scatter_mean_',...
        num2str(nruns), 'runs'));
    save_figure(h, 9, 8, fname, '.pdf');
end

% old version
%{
h = figure(5);
hold on
grid on
grid minor
scatter(x(:), y(:)/4, 100);
%plot(K12_all, periods_stat(:, 1), 'rd', 'MarkerSize', 20);
%errorbar(N, periods_stat(:, 1), periods_stat(:, 2), 'ro');
%set(gca, 'YScale', 'log');
xlim([K12_all(1)-1 K12_all(end)+1]);
xlabel('$$K_{12}$$');
ylabel('Period');
set(gca,'FontSize', 24);
%set(gca,'YScale','log');
set(h, 'Units', 'Inches', 'Position', [1 1 9 8]);
%}
%% Histograms (not used)
%% Plot t_out histogram
% (for single N value at a time)
%for idxN=1:nN
    %idxN = 1; % 
    %thisN = N(idxN);
    %t_data = t_out(idxN, :, :);

% Plot together
h = figure(43);
hold on
%tmax = ceil(max(t_data(:))/10)*10;
%edges = 0:tmax/20:tmax;
nbins = 10;
for i=1:numel(K_all)
    histogram(t_out(i,:), nbins, 'Normalization', 'probability');
end
xlabel('$$t_{out}$$');
ylabel('probability');
%xlim([0 tmax]);
set(gca,'FontSize', 24);
set(gca, 'XScale', 'log');
set(h, 'Units', 'Inches', 'Position', [1 1 9 8]);
legend(sprintfc('$$K_{12}=%d$$',K_all), 'Interpreter', 'latex');

qsave = 1;
if qsave
    fname = fullfile(save_path_fig, strcat(fname_str, '_t_out_hist_all_',...
        num2str(nruns) ,'runs'));
    save_figure(h, 9, 8, fname, '.pdf');
end
    %close all
%end
%% Plot period histogram
% (for single N value at a time)
close all
for idxK=1:nK
    %idxN = 5; % 
    %thisN = N(idxN);
    K22 = K_all(idxK);
    p_data = periods(idxK, :);
    uniq = unique(p_data);
    C = categorical(p_data, uniq); 
    h1 = figure(idxK);
    histogram(C, 'Normalization', 'count');
    xlabel('period');
    ylabel('count');
    title(sprintf('K12=%d, %d simulations', K22, nruns));
    set(gca,'FontSize', 24);
    set(h1, 'Units', 'Inches', 'Position', [1 1 9 8]);

    qsave = 1;
    save_path = save_path_fig;
    %save_path = 'H:\My Documents\Multicellular automaton\latex\10_two_signals\analyze_trajectories\analyze_trajectories_vs_K12\all_period_histograms\randomized_lattice_N225';
    
    if qsave
        fname = fullfile(save_path, strcat('K_',...
            num2str(K_all(idxK)), '_period_hist'));
        save_figure(h1, 9, 8, fname, '.pdf');
    end
    %close all
end
