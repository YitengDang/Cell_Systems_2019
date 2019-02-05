clear all
close all
clc
set(0, 'defaulttextinterpreter', 'latex')
%% Parameters
%parent_folder = 'H:\My Documents\Multicellular automaton\data\two_signals\parameter_gradient';
%parent_folder = 'L:\BN\HY\Shared\Yiteng\two_signals\parameter_gradient_sinusoidal';
parent_folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\parameter_gradient_step_function';
%subfolder = 'negative_control';
%subfolder = 'vertical_step_function_Ay_0p5_ny_1';

wave_type = 'square_wave';
%wave_type = 'sine_wave';

Ax_all = 0.1:0.1:0.5;
Ay_all = 0.1:0.1:0.5;

Ax = 0;
nx = 1;
Ay = 0;
ny = 1;
%% Process raw data

orientation = 'Ax';
if strcmp(orientation, 'Ax')
    A_all = Ax_all;
elseif strcmp(orientation, 'Ay')
    A_all = Ay_all;
end
for folder_count=5 %1:numel(A_all)
    if strcmp(orientation, 'Ax')
        Ax = A_all(folder_count);
        fprintf('Ax = %.1f \n', Ax);
        if strcmp(wave_type, 'square_wave')
            subfolder = strrep(...
                sprintf('horizontal_step_function_Ax_%.1f_nx_%d', Ax, nx),...
                '.', 'p');
            pattern = 'Parameter_gradient_K_2_1_square_wave_t_out_(\d+)_period_(Inf|\d+)-v\d+'; % for regexp to recognize files
        end
    elseif strcmp(orientation, 'Ay')
        Ay = A_all(folder_count);
        fprintf('Ay = %.1f \n', Ay);
        if strcmp(wave_type, 'square_wave')
            subfolder = strrep(...
                sprintf('vertical_step_function_Ay_%.1f_ny_%d', Ay, nx),...
                '.', 'p');
            pattern = 'Parameter_gradient_K_2_1_square_wave_t_out_(\d+)_period_(Inf|\d+)-v\d+'; % for regexp to recognize files
        end
    end
   
    if strcmp(wave_type, 'sine_wave')
        subfolder = strrep(...
            sprintf('sine_wave_gradient_K21_Ax_%.1f_nx_%d_Ay_%.1f_ny_%d',...
            Ax, nx, Ay, ny), '.', 'p');
        pattern = 'Parameter_gradient_K_2_1_square_wave_t_out_(\d+)_period_(Inf|\d+)-v\d+';  % for regexp to recognize files
    end
    
    folder = fullfile(parent_folder, subfolder);
    listing = dir(folder);
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

    %%
    orientation_all = {};
    bands_in_wave_all = [];
    number_of_waves_all = [];
    diag_vec_all = {};
    band_vec_all = {};
    wave_state_all = [];
    bended_all = [];
    true_trav_wave = []; % records whether a found wave is "true" travelling wave as determined by travelling_wave_test
    %trav_wave = [];
    nruns = 0; % total number of simulations analyzed
    fname_all = {};
    
    for file_count = 1:num_files
        disp(file_count);
        nruns = nruns+1;

        % skip non-periodic files
        filename = names{file_count};
        [tokens, out] = regexp( filename, pattern, 'tokens', 'match');
        if str2double(tokens{1}{2})==Inf
            disp('Period = Inf, continue');
            continue
        end

        % load file
        fname = fullfile(folder, names{file_count});
        load(fname);
        a0 = save_consts_struct.a0;

        % Find periodicity and t_onset
        [period_ub, ~] = periodicity_test_short(cells_hist);
        [period, t_onset] = periodicity_test_detailed(cells_hist, t_out,...
            period_ub);
        t_wave = t_onset; 

        % Find whether it is a travelling wave
        [trav_wave, trav_wave_2]  = travelling_wave_test(cells_hist, a0,...
            period, t_out, distances);
        if trav_wave
            % Find wave properties
            wave_state = -1;
            [orientation, bands_in_wave, number_of_waves, diag_vec, band_vec,...
                wave_state, bended] = determine_wave_properties(cells_hist, t_wave, wave_state);
            %[wave_state_translated] = translate_states(wave_state);
        else
            % otherwise, the condition is that the predictor should give the
            % same direction for all configurations during one period
            orientations_this_sim = {};
            try 
                wave_state = -1;
                [orientation, bands_in_wave, number_of_waves, diag_vec, band_vec,...
                    wave_state, bended] = determine_wave_properties(cells_hist, t_wave, wave_state);
                [wave_state_translated] = translate_states(wave_state);

                orientations_this_sim{end+1} = orientation;
            catch
                disp('Not a travelling wave');
                continue
            end
            if numel(unique(orientations_this_sim))>1
                disp('Multiple orientations, continue');
                continue            
            end
        end
        
        % Store properties
        orientation_all{end+1} = orientation;
        bands_in_wave_all(end+1) = bands_in_wave;
        number_of_waves_all(end+1) = number_of_waves;
        diag_vec_all{end+1} = diag_vec;
        band_vec_all{end+1} = band_vec;
        wave_state_all(end+1) = wave_state;
        bended_all(end+1) = bended;
        true_trav_wave(end+1) = trav_wave;
        %trav_wave(end+1) = trav_wave_2;
        fname_all{end+1} = fname;
        %}
    end
    
    %% Save analyzed data
    %
    fname_str = sprintf('analyzed_data_%s_%druns', subfolder, nruns);
    fname = fullfile(parent_folder, fname_str);
    save(fname, 'subfolder', 'orientation_all', 'bands_in_wave_all',...
      'number_of_waves_all', 'diag_vec_all', 'band_vec_all',...
      'wave_state_all', 'bended_all', 'true_trav_wave',...
      'nruns', 'fname_all');
    %}
end

%% Save analyzed data
%{
fname_str = sprintf('analyzed_data_%s_%druns', subfolder, nruns);
fname = fullfile(parent_folder, fname_str);
save(fname, 'subfolder', 'orientation_all', 'bands_in_wave_all',...
    'number_of_waves_all', 'diag_vec_all', 'band_vec_all',...
    'wave_state_all', 'bended_all', 'true_trav_wave', 'nruns');
%}
%% Load specific saved data
%
%parent_folder = 'H:\My Documents\Multicellular automaton\data\two_signals\parameter_gradient';
parent_folder= 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\parameter_gradient_step_function';
%subfolder = 'negative_control';
Ax = 0.1;
nx = 1;
subfolder = strrep(...
        sprintf('horizontal_step_function_Ax_%.1f_nx_%d', Ax, nx),...
        '.', 'p');
    
nruns = 200;
fname_str = sprintf('analyzed_data_%s_%druns', subfolder, nruns);
fname = fullfile(parent_folder, fname_str);
load(fname);
%}

% find the names of the trajectories generating travelling waves;
%[idx] = find(true_trav_wave);
%fname_all(idx(1))
%fname_all(idx(2))
%fname_all(idx(3))
%}

% copy "true" and "false" TWs to separate folder
fname_strict_TW = fname_all(~~true_trav_wave);
fname_loose_TW = fname_all(~true_trav_wave);
folder_strict = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\parameter_gradient_step_function\TW_strict_filtered';
folder_loose = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\parameter_gradient_step_function\TW_loose_filtered';

for ii=1:numel(fname_strict_TW)
    fname = strcat(fname_strict_TW{ii}, '.mat');
    [~, name, ext] = fileparts(fname);
    fname_str = strrep(sprintf('Ax_%.1f_nx_%d_%s', Ax, nx, name), '.', 'p');
    fname_new = fullfile(folder_strict, strcat(fname_str, ext) );
    copyfile(fname, fname_new);
end

for ii=1:numel(fname_loose_TW)
    fname = strcat(fname_loose_TW{ii}, '.mat');
    [~, name, ext] = fileparts(fname);
    fname_str = strrep(sprintf('Ax_%.1f_nx_%d_%s', Ax, nx, name), '.', 'p');
    fname_new = fullfile(folder_loose,  strcat(fname_str, ext) );
    copyfile(fname, fname_new);
end
%% Analyze found waves
%{
% Orientation
C = categorical(orientation_all);

h = figure;
histogram(C)
ylim([0 100]);

qsave = 0;
folder = 'H:\My Documents\Multicellular automaton\figures\parameter_gradient';
fname_str = sprintf('orientation_histogram_%s', subfolder);
fname = fullfile(folder, fname_str);
save_figure(h, 10, 8, fname, '.pdf', qsave)
%}

%% Plot multiple ones at the same time
% -----(1) for square waves --------------
%
% folder with analyzed data
parent_folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\parameter_gradient_step_function'; 
nruns = 200;
Avals = 0.1:0.1:0.5;

%
% --- vertical ---
ny = 1;
subfolders_temp = strrep(strcat(...
    sprintfc('vertical_step_function_Ay_%.1f_', Avals), ...
    'ny_', num2str(ny)), '.', 'p');
orientations = {'Horizontal', 'Vertical', 'Diagonal'};
grad = 'Ay';
%}
% --- horizontal ---
%{
nx = 1;
subfolders_temp = strrep(strcat(...
    sprintfc('horizontal_step_function_Ax_%.1f_', Avals), ...
    'nx_', num2str(nx)), '.', 'p');
orientations = {'Vertical', 'Horizontal','Diagonal'};
grad = 'Ax';
%}
% ------------------
%subfolders_temp = strrep(sprintfc('vertical_step_function_Ay_%.1f_ny_%d',...
%    [0.1:0.1:0.5]), '.', 'p');
subfolders = [{'negative_control'}, subfolders_temp];
y_counts = zeros(numel(subfolders), 3);
y_counts_true = zeros(numel(subfolders), 3);
for i=1:numel(subfolders)
    disp(i);
    subfolder = subfolders{i};
    fname_str = sprintf('analyzed_data_%s_%druns', subfolder, nruns);
    fname = fullfile(parent_folder, fname_str);
    load(fname, 'orientation_all', 'true_trav_wave');
    disp(numel(orientation_all));
    
    % Do for "all" travelling waves including those that do not meet the
    % strict criteria of travelling_wave_test
    y_counts(i, 1) = sum(strcmp(orientation_all, orientations{1}))/nruns; 
    y_counts(i, 2) = sum(strcmp(orientation_all, orientations{2}))/nruns; 
    y_counts(i, 3) = sum(strcmp(orientation_all, orientations{3}))/nruns; 
    
    % Also do for "true travelling waves"
    idx = find(true_trav_wave);
    if ~isempty(idx)
        orientation_all_true = orientation_all(idx);
        y_counts_true(i, 1) = sum(strcmp(orientation_all_true, orientations{1}))/nruns; 
        y_counts_true(i, 2) = sum(strcmp(orientation_all_true, orientations{2}))/nruns; 
        y_counts_true(i, 3) = sum(strcmp(orientation_all_true, orientations{3}))/nruns; 
    end
end
%}
%%
% -----(2) for sine waves --------------
%{
parent_folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\parameter_gradient_step_function'; %'L:\BN\HY\Shared\Yiteng\two_signals\parameter_gradient_sinusoidal';
nruns = 100;
Avals = 0.1:0.1:0.5;

orientations = {'Vertical', 'Horizontal','Diagonal'};
grad = 'Ax';

if strcmp(gradient, 'Ax')
    % --- (Ax gradient) ---
    %subfolders = strrep(...
    %    sprintfc('sine_wave_gradient_K21_Ax_%.1f_nx_1_Ay_0p0_ny_1', Avals), '.', 'p');
elseif strcmp(gradient, 'Ay')
    % --- (Ay gradient) ---
    subfolders = strrep(...
       sprintfc('sine_wave_gradient_K21_Ax_0p0_nx_1_Ay_%.1f_ny_1',...
       Avals), '.', 'p');
    % ------------------
end

y_counts = zeros(numel(subfolders), 3);
y_counts_true = zeros(numel(subfolders), 3);
for i=1:numel(subfolders)
    disp(i);
    subfolder = subfolders{i};
    fname_str = sprintf('analyzed_data_%s_%druns', subfolder, nruns);
    fname = fullfile(parent_folder, fname_str);
    load(fname, 'orientation_all', 'true_trav_wave'); %, 'num_simulations');
    disp(numel(orientation_all));
    
    % Do for "all" travelling waves including those that do not meet the
    % strict criteria of travelling_wave_test
    y_counts(i, 1) = sum(strcmp(orientation_all, orientations{1}))/nruns;
    y_counts(i, 2) = sum(strcmp(orientation_all, orientations{2}))/nruns;
    y_counts(i, 3) = sum(strcmp(orientation_all, orientations{3}))/nruns;
    
    % Also do for "true travelling waves"
    idx = find(true_trav_wave);
    if ~isempty(idx)
        orientation_all_true = orientation_all(idx);
        y_counts_true(i, 1) = sum(strcmp(orientation_all_true, orientations{1}))/nruns; 
        y_counts_true(i, 2) = sum(strcmp(orientation_all_true, orientations{2}))/nruns; 
        y_counts_true(i, 3) = sum(strcmp(orientation_all_true, orientations{3}))/nruns; 
    end
end
%}
%%
h1 = figure;
b = bar(y_counts, 'stacked', 'FaceColor', 'flat');
legend(orientations);
if strcmp(grad, 'Ax')
    %xlabel('Relative strength of x gradient ($$A_x$$)');
    new_order = [2 1 3]; % new order of default colors
elseif strcmp(grad, 'Ay')
    xlabel('Relative strength of y gradient ($$A_y$$)');
    new_order = [1 2 3]; % new order of default colors
end
ylabel('Frequency');
title('All periodic solutions');
%set(gca,'xticklabel', [{'0 (control)'}, sprintfc('%.1f', Avals)]);
set(gca,'xticklabel', sprintfc('%.1f', Avals));
set(gca, 'FontSize', 18);
ylim([0 1]);
% set colors of bar
cmp = get(gca,'colororder');
for k = 1:3
    b(k).CData = cmp(new_order(k), :);
end

h2 = figure;
b=bar(y_counts_true, 'stacked', 'FaceColor', 'flat');
legend(orientations);
if strcmp(grad, 'Ax')
    %xlabel('Relative strength of x gradient ($$A_x$$)');
    new_order = [2 1 3]; % new order of default colors
elseif strcmp(grad, 'Ay')
    xlabel('Relative strength of y gradient ($$A_y$$)');
    new_order = [1 2 3]; % new order of default colors
end
ylabel('Frequency');
title('(Real) travelling waves');
%set(gca,'xticklabel', [{'0 (control)'}, sprintfc('%.1f', Avals)]);
set(gca,'xticklabel', sprintfc('%.1f', Avals));
set(gca, 'FontSize', 18);
ylim([0 1]);
% set colors of bar
cmp = get(gca, 'colororder');
%new_order = [1 2 3]; % new order of default colors
for k = 1:3
    b(k).CData = cmp(new_order(k), :);
end

% save figures
qsave = 0;
folder = 'H:\My Documents\Multicellular automaton\figures\parameter_gradient';
subfolder = 'step_gradient';
if strcmp(grad, 'Ax')
    fname_str1 = sprintf('orientation_vs_gradient_x_bar_stacked_all_waves_nx_%d_%druns', nx, nruns);
    fname_str2 = sprintf('orientation_vs_gradient_x_bar_stacked_true_trav_waves_nx_%d_%druns', nx, nruns);
elseif strcmp(grad, 'Ay')
    fname_str1 = sprintf('orientation_vs_gradient_y_bar_stacked_all_waves_ny_%d_%druns', ny, nruns);
    fname_str2 = sprintf('orientation_vs_gradient_y_bar_stacked_true_trav_waves_ny_%d_%druns', ny, nruns);
end
save_figure(h1, 10, 8, fullfile(folder, subfolder, fname_str1), '.pdf', qsave)
save_figure(h2, 10, 8, fullfile(folder, subfolder, fname_str2), '.pdf', qsave)
%}
%% Number of waves
%{
h = figure;
histogram(number_of_waves_all);
ylim([0 numel(number_of_waves_all)]);
xlabel('Number of waves');

%% Bands in waves
h = figure;
histogram(bands_in_wave_all);
xlabel('Number of bands');

%% Bended waves
h = figure;
C = categorical(bended_all, [0 1], {'Straight', 'Bent'});
histogram(C);
%}

%%
