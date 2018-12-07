% Analyze batch simulation results for moving cells
close all
clear all
% maxNumCompThreads(4);
% warning off
set(0, 'defaulttextinterpreter', 'latex');

%% Set folders
%parent_folder = 'W:\staff-bulk\tnw\BN\HY\Shared\Yiteng\two_signals\moving_cells';
parent_folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\moving_cells';

% Save data folder
save_folder = parent_folder;

% Save figures folder
save_fig_folder = parent_folder;
%% One signal
%
sigma_D_all = 0.01; %[0.001 0.01 0.02 0.05 0.1];
nruns = 20;
t_max = 10000;
gz = 15;
N = gz^2;
num_cells_reached_all = cell(numel(sigma_D_all), nruns); 
p_final_all = zeros(numel(sigma_D_all), nruns);

subfolder = 'one_signal_temp';
folder = fullfile(parent_folder, subfolder);

for ii=1:numel(sigma_D_all)
    sigma_D = sigma_D_all(ii);
    for jj=1:nruns
        %fname_str = 'one_signal_sigma_D_0p001_t_out_1000-v1';
        fname_str = strrep(sprintf('one_signal_sigma_D_%.3f_t_out_%d-v%d',...
            sigma_D, t_max, jj), '.' ,'p');
        fname = fullfile(folder, strcat(fname_str, '.mat'));
        if exist(fname, 'file')==2
            disp(fname);
            load(fname)
            
            cells_final = cells_hist{end};
            p_final_all(ii,jj) = sum(cells_final)/N;
            
            %
            % number of cells where signal has reached over time
            num_cells_reached = zeros(numel(cells_hist), 1);
            
            % indices of cells where the signal has reached
            cells_reached = sum(cells_hist{1}, 2)>0; 
            num_cells_reached(1) = sum(cells_reached);
            for tt=2:numel(cells_hist)
                reached_temp = sum(cells_hist{tt}, 2)>0;
                % update newly reached cells
                cells_reached(reached_temp) = 1; 
                % determine number of cells reached
                num_cells_reached(tt) = sum(cells_reached);
            end
            num_cells_reached_all{ii,jj} = num_cells_reached;
        end
    end
end

%% Save analyzed data
save_fname_str = sprintf('analyzed_data_%s_t_max_%d_nruns_%d',...
    subfolder, t_max, nruns);
save_file = fullfile(parent_folder, save_fname_str);
save(save_file, 'subfolder', 'sigma_D_all', 'p_final_all', 'num_cells_reached_all');

%% Display all p_final
%{
h=figure;
x_data = repmat(sigma_D_all', 1, nruns);
scatter(x_data(:), p_final_all(:))
set(gca, 'XScale', 'log');
%}

%% Plot fractions over time

% Calculate avg. trajectory
counter = zeros( numel(sigma_D_all), 1 );
avg_trajectories_all = zeros( numel(sigma_D_all), t_max+1);
for ii=1:numel(sigma_D_all)
    for jj=1:nruns
        if ~isempty(num_cells_reached_all{ii,jj})
            avg_trajectories_all(ii, :) = avg_trajectories_all(ii, :) +...
                (num_cells_reached_all{ii,jj}/N)';
            counter(ii) = counter(ii) + 1;
        end
    end
end
avg_trajectories_all = avg_trajectories_all./counter;

% plot 
colors = {'b', 'r', 'c', 'm', 'g'};
h = figure;
p_temp = cell(size(colors));
hold on
% Calculate fractions and plot
for ii=1:numel(sigma_D_all)
    p_temp{ii} = plot(Inf, Inf, colors{ii});
    for jj=1:size(num_cells_reached_all, 2)
        if ~isempty(num_cells_reached_all{ii,jj})
            frac_cells_reached = num_cells_reached_all{ii,jj}/N;
            plot(0:t_max, frac_cells_reached, colors{ii});
        end
    end
    
    % plot avg. trajectory
    plot(0:t_max, avg_trajectories_all(ii, :), colors{ii}, 'LineWidth', 3);
end
legend([p_temp{:}], sprintfc('$\\sigma_D = %.3f$', sigma_D_all),...
    'Interpreter', 'latex', 'Location', 'eo');
ylim([0.5 1]);
xlabel('Time');
ylabel('Fraction ON cells');
set(gca, 'FontSize', 24);

% Save figure
qsave = 1;
% save_folder = 'W:\staff-bulk\tnw\BN\HY\Shared\Yiteng\two_signals\moving_cells';
fname_str = sprintf('analyzed_data_%s_t_max_%d_nruns_%d_p_vs_t',...
    subfolder, t_max, nruns);
fname = fullfile(save_fig_folder, fname_str);
save_figure(h, 10, 8, fname, '.pdf', qsave);

% plot average trajectory only
%{
figure;
plot(1:t_max+1, avg_trajectories_all);
%}
%% Two signals
% 
% Subdomain oscillations -> Determine how many of the cells the signal has reached
% Period 4 oscillations -> Look at I(t) or Theta(t)
%sigma_D_all = [0.001 0.01 0.1];
%sigma_D_all = [0.01 0.1 1];
sigma_D_all = 0.001;
nruns = 20;
t_max = 2000;

%parent_folder = 'W:\staff-bulk\tnw\BN\HY\Shared\Yiteng\two_signals\moving_cells';
subfolder = 'subdomain_oscillations';
%subfolder = 'period_4_oscillations';

% load test file to get parameters
sigma_D = 0.001; jj = 1;
%fname_str = strrep(sprintf('two_signal_mult_sigma_D_%.3f_t_out_%d_period_Inf-v%d',...
%	sigma_D, t_max, jj), '.' ,'p');
fname_str = strrep(sprintf('two_signal_mult_sigma_D_%.3f_t_out_%d-v%d',...
	sigma_D, t_max, jj), '.' ,'p');
fname = fullfile(parent_folder, subfolder, strcat(fname_str, '.mat'));
load(fname, 'save_consts_struct')
a0 = save_consts_struct.a0;
N = save_consts_struct.N;
gz = sqrt(N);

%% Calculate from data 
% Variables to store
% number of cells that have turned at least one gene on (assuming background = (0,0))
% num_cells_reached_all = cell(numel(sigma_D_all), nruns); 
theta_all = cell(numel(sigma_D_all), nruns); 
I_all = theta_all;
% Load all data
for ii=1:numel(sigma_D_all)
    sigma_D = sigma_D_all(ii);
    for jj=1:nruns
        %fname_str = 'one_signal_sigma_D_0p001_t_out_1000-v1';
        %fname_str = strrep(sprintf('two_signal_mult_sigma_D_%.3f_t_out_%d_period_Inf-v%d',...
        %    sigma_D, t_max, jj), '.' ,'p');
        fname_str = strrep(sprintf('two_signal_mult_sigma_D_%.3f_t_out_%d-v%d',...
            sigma_D, t_max, jj), '.' ,'p');
        
        fname = fullfile(parent_folder, subfolder, strcat(fname_str, '.mat'));
        if exist(fname, 'file')==2
            disp(fname);
            load(fname)
            %------------------------------------------------------------
            % period 4 oscillations
            %{
            theta_t = zeros( numel(cells_hist), 2 );
            I_t = zeros( numel(cells_hist), 2 );
            for tt=1:numel(cells_hist)
                pos = positions_all{tt};
                Lx = 1;
                delx = Lx/gz;
                dely = sqrt(3)/2*delx;
                Ly = dely*gz;
                dist = calc_dist_periodic(pos(:,1), pos(:,2), Lx, Ly);
                [I_t(tt,1), theta_t(tt, 1)] = moranI(cells_hist{tt}(:,1), a0*dist);
                [I_t(tt,2), theta_t(tt, 2)] = moranI(cells_hist{tt}(:,2), a0*dist);
            end
            theta_all{ii, jj} = theta_t;
            I_all{ii, jj} = I_t;
            %}
            %------------------------------------------------------------
            % subdomain oscillations
            %
            % number of cells where signal has reached over time
            num_cells_reached = zeros(numel(cells_hist), 1);
            
            % indices of cells where the signal has reached
            cells_reached = sum(cells_hist{1}, 2)>0; 
            num_cells_reached(1) = sum(cells_reached);
            for tt=2:numel(cells_hist)
                reached_temp = sum(cells_hist{tt}, 2)>0;
                % update newly reached cells
                cells_reached(reached_temp) = 1; 
                % determine number of cells reached
                num_cells_reached(tt) = sum(cells_reached);
            end
            num_cells_reached_all{ii,jj} = num_cells_reached;
            %}
        end
    end
end

%% Save analyzed data

%save_folder = parent_folder;
save_fname_str = sprintf('analyzed_data_%s_t_max_%d_nruns_%d_sigma_D_0p001',...
    subfolder, t_max, nruns);
save_file = fullfile(save_folder, save_fname_str);
%{
% Save period 4 osc. data
% save(save_file, 'subfolder', 'sigma_D_all', 'theta_all', 'I_all');

% Save subdomain osc. data
save(save_file, 'subfolder', 'sigma_D_all', 'num_cells_reached_all');
%}
%% Load analyzed data
save_fname_str = sprintf('analyzed_data_%s_t_max_%d_nruns_%d_incomplete',...
    subfolder, t_max, nruns);

% Load subdomain oscillations
%{
fname_str = 'analyzed_data_subdomain_oscillations_t_max_1000_nruns_20';
save_file = fullfile(parent_folder, fname_str);
load(save_file, 'subfolder', 'sigma_D_all', 'num_cells_reached_all');
%}
% Load period 4 oscillations
fname_str = 'analyzed_data_period_4_oscillations_t_max_20_nruns_incomplete';
save_file = fullfile(save_folder, fname_str);
load(save_file, 'subfolder', 'sigma_D_all', 'theta_all', 'I_all');
%}
%% Plot fractions over time
% register fractions at end of simulation
frac_end = zeros(size(num_cells_reached_all));

% plot 
colors = {'b', 'r', 'g'};
h = figure;
p_temp = {};
hold on
for ii=1 %:numel(sigma_D_all)
    p_temp{ii} = plot(Inf, Inf, colors{ii});
    for jj=1:size(num_cells_reached_all, 2)
        frac_cells_reached = num_cells_reached_all{ii,jj}/N;
        plot(1:numel(frac_cells_reached), frac_cells_reached, colors{ii});
        frac_end(ii,jj) = frac_cells_reached(end);
    end
end
legend([p_temp{:}], sprintfc('$\\sigma_D = %.3f$', sigma_D_all),...
    'Interpreter', 'latex', 'Location', 'se');

ylim([0.4 1]);
xlabel('Time');
ylabel('Fraction cells reached');
set(gca, 'FontSize', 24);

% Save figure
qsave = 1;
% save_folder = 'W:\staff-bulk\tnw\BN\HY\Shared\Yiteng\two_signals\moving_cells';
fname_str = sprintf('analyzed_data_%s_t_max_%d_nruns_%d_frac_reached_vs_t_sigma_D_0p001',...
    subfolder, t_max, nruns);
fname = fullfile(save_fig_folder, fname_str);
save_figure(h, 10, 8, fname, '.pdf', qsave);
%}
% Plot fractions at end
%{
h=figure;
hold on
bins = 0:0.1:1;
for ii=1:numel(sigma_D_all)
    histogram(frac_end(ii,:), bins )
end
%}
%% Plot Theta(t) for different trajectories
colors = {'b', 'r', 'g'};
h1 = figure;
p_temp = {};
hold on
for ii=1:numel(sigma_D_all)
    p_temp{ii} = plot(Inf, Inf, 'x', 'Color', colors{ii});
    for jj=1:size(theta_all, 2)
        if ~isempty(theta_all{ii,jj})
            plot(0:t_max, theta_all{ii,jj}(:,1),...
                'x', 'Color', colors{ii});
        end
    end
end
legend([p_temp{:}], sprintfc('$\\sigma_D = %.3f$', sigma_D_all),...
    'Interpreter', 'latex', 'Location', 'eo');
xlabel('$t$');
ylabel('$\Theta^{(1)}(t)/f_N$');
xlim([0 t_max]);
ylim([-0.2 1]);
set(gca, 'FontSize', 24);

%----
colors = {'b', 'r', 'g'};
h2 = figure;
p_temp = {};
hold on
for ii=1:numel(sigma_D_all)
    p_temp{ii} = plot(Inf, Inf, 'x', 'Color', colors{ii});
    for jj=1:size(theta_all, 2)
        if ~isempty(theta_all{ii,jj})
            plot(0:t_max, theta_all{ii,jj}(:, 2),...
                'x', 'Color', colors{ii});
        end
    end
end
legend([p_temp{:}], sprintfc('$\\sigma_D = %.3f$', sigma_D_all),...
    'Interpreter', 'latex', 'Location', 'eo');
xlabel('$t$');
ylabel('$\Theta^{(2)}(t)/f_N$');
xlim([0 t_max]);
ylim([-0.2 1]);
set(gca, 'FontSize', 24);

% Save figures
qsave = 1;
% save_folder = 'W:\staff-bulk\tnw\BN\HY\Shared\Yiteng\two_signals\moving_cells';
fname_str = sprintf('analyzed_data_%s_t_max_%d_nruns_%d_theta1_vs_t',...
    subfolder, t_max, nruns);
fname = fullfile(save_fig_folder, fname_str);
save_figure(h1, 10, 8, fname, '.pdf', qsave);

fname_str = sprintf('analyzed_data_%s_t_max_%d_nruns_%d_theta2_vs_t',...
    subfolder, t_max, nruns);
fname = fullfile(save_fig_folder, fname_str);
save_figure(h2, 10, 8, fname, '.pdf', qsave);
%%
%{
figure;
ii = 1;
jj = 1;
plot(1:size(theta_all{ii,jj}, 1), theta_all{ii,jj}(:,1), 'o', 'Color', colors{ii});
%}