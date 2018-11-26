%% Analyze trajectory output times and periods
clear all
close all
set(0, 'defaulttextinterpreter', 'latex');
%%
path = 'L:\BN\HY\Shared\Yiteng\two_signals\parameter set 2b';
%path = 'L:\BN\HY\Shared\Yiteng\two_signals\sweep K12';
%N = [5 8 10 12 13 15 20].^2;
gz = 15;
N = gz^2;
nN = numel(N);

Ilist = 0:0.1:0.7; 
%{[0 0], [0.1 0.1], [0.2 0.2], [0.3 0.3], [0.4 0.4], [0.5 0.5], [0.6 0.6], [0.7 0.7]};

folders = sprintfc('N%d',N);
%folder = 'N225'; % 'N64', 'N100', 'N225';
%fname_str = sprintf('parameter_set_2b_%s', folder);
fname_str = 'parameter_set_2b';
save_path_fig = 'H:\My Documents\Multicellular automaton\figures\two_signals\analyze_trajectories_final_pI_stats';

a0 = 1.5;
nruns = 100;
I_ini = zeros(nN, numel(Ilist), nruns, 2);
I1_out = zeros(nN, numel(Ilist), nruns);
I2_out = zeros(nN, numel(Ilist), nruns);
p1_out = zeros(nN, numel(Ilist), nruns);
p2_out = zeros(nN, numel(Ilist), nruns);
periods = zeros(nN, numel(Ilist), nruns);

%names = {};
for idxN=1:nN
    gz = sqrt(N(idxN));
    dist = init_dist_hex(gz, gz);
    
    folder = folders{idxN};
    listing = dir(fullfile(path, folder));
    num_files = numel(listing)-2; %first two entries are not useful

    d = '(\d+)';
    d1 = '0p(\d)0';
    d2 = '(\d+|Inf)';
    %pattern = sprintf('two_signal_mult_%s_chaotic_state_search_t_out_%s_period_%s-v%s',...
    %    folder, d, d2, d); 
    for idxI=1:numel(Ilist)
        I_ini_temp = Ilist(idxI)*ones(2,1);
        pattern = strrep(sprintf('two_signal_mult_%s_initiateI1_I_ini_%.2f_%.2f_t_out_%s_period_%s-v%s',...
            folder, I_ini_temp(1), I_ini_temp(2), d, d2, d), '.', 'p');
        %pattern = strrep(sprintf('two_signal_mult_N%d_initiateI0_K12_%d_t_out_%s_period_%s-v%s',...
        %    N, K12, d, d2, d), '.', 'p');
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
                    cells_out = cells_hist{end};
                    I_ini(idxN, idxI, count+1, 1) = moranI(cells_ini(:, 1), a0*dist);
                    I_ini(idxN, idxI, count+1, 2) = moranI(cells_ini(:, 2), a0*dist);                
                    p1_out(idxN, idxI, count+1) = mean(cells_out(:,1), 1);
                    p2_out(idxN, idxI, count+1) = mean(cells_out(:,2), 1);
                    I1_out(idxN, idxI, count+1) = moranI(cells_out(:, 1), a0*dist);
                    I2_out(idxN, idxI, count+1) = moranI(cells_out(:, 2), a0*dist);
                    periods(idxN, idxI, count+1) = str2double(tokens{1}{2});
                    count = count + 1;
                    %names{count} = name;
                end
            end
        end

    end

end

%% Load analyzed data


%% Plot statistics
% p1-p2
xdata = p1_out(periods~=Inf);
ydata = p2_out(periods~=Inf);

h1 = figure(1);
hold on
scatter(xdata, ydata);
xlabel('final $$p_1$$');
ylabel('final $$p_2$$');
set(gca,'FontSize', 24);
set(h1, 'Units', 'Inches', 'Position', [1 1 9 8]);
legend(folders, 'Location', 'se');
xlim([0 1]);
ylim([0 1]);

qsave = 0;
if qsave
    fname = fullfile(save_path_fig, strcat(fname_str, '_'));
    save_figure(h1, 9, 8, fname, '.pdf');
end

% p1-I1
xdata = p1_out(periods~=Inf);
ydata = I1_out(periods~=Inf);

h2 = figure(2);
hold on
scatter(xdata, ydata);
xlabel('final $$p_1$$');
ylabel('final $$I_1$$');
set(gca,'FontSize', 24);
set(h2, 'Units', 'Inches', 'Position', [1 1 9 8]);
legend(folders, 'Location', 'se');
xlim([0 1]);
ylim([-1 1]);

qsave = 0;
if qsave
    fname = fullfile(save_path_fig, strcat(fname_str, '_'));
    save_figure(h2, 9, 8, fname, '.pdf');
end

% p2-I2
xdata = p2_out(periods~=Inf);
ydata = I2_out(periods~=Inf);

h3 = figure(3);
hold on
scatter(xdata, ydata);
xlabel('final $$p_2$$');
ylabel('final $$I_2$$');
set(gca,'FontSize', 24);
set(h3, 'Units', 'Inches', 'Position', [1 1 9 8]);
legend(folders, 'Location', 'se');
xlim([0 1]);
ylim([-1 1]);

qsave = 0;
if qsave
    fname = fullfile(save_path_fig, strcat(fname_str, '_'));
    save_figure(h3, 9, 8, fname, '.pdf');
end

% I1-I2
xdata = I1_out(periods~=Inf);
ydata = I2_out(periods~=Inf);

h4 = figure(4);
hold on
scatter(xdata, ydata);
xlabel('final $$I_1$$');
ylabel('final $$I_2$$');
set(gca,'FontSize', 24);
set(h4, 'Units', 'Inches', 'Position', [1 1 9 8]);
legend(folders, 'Location', 'se');
xlim([-1 1]);
ylim([-1 1]);

qsave = 0;
if qsave
    fname = fullfile(save_path_fig, strcat(fname_str, '_'));
    save_figure(h4, 9, 8, fname, '.pdf');
end