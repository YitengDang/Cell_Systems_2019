% Analyze batch simulation results for moving cells
close all
clear all
% maxNumCompThreads(4);
% warning off
%% One signal
%
sigma_D_all = [0.001 0.01 0.1];
nruns = 20;
t_max = 1000;
N = 15^2;
num_cells_reached_all = cell(numel(sigma_D_all), nruns); 
p_final_all = zeros(numel(sigma_D_all), nruns);
folder = 'W:\staff-bulk\tnw\BN\HY\Shared\Yiteng\two_signals\moving_cells\one_signal';

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

%% Display all p_final
h=figure;
x_data = repmat(sigma_D_all', 1, nruns);
scatter(x_data(:), p_final_all(:))
set(gca, 'XScale', 'log');
%% Plot fractions over time
colors = {'b', 'r', 'g'};
h = figure;
p_temp = {};
hold on
for ii=1:numel(sigma_D_all)
    p_temp{ii} = plot(Inf, Inf, colors{ii});
    for jj=1:size(num_cells_reached_all, 2)
        frac_cells_reached = num_cells_reached_all{ii,jj}/N;
        plot(1:numel(frac_cells_reached), frac_cells_reached, colors{ii});
    end
end
legend([p_temp{1} p_temp{2} p_temp{3}], sprintfc('$\\sigma_D = %.3f$', sigma_D_all),...
    'Interpreter', 'latex');
%}
%% Two signals
%{
% Subdomain oscillations
% Determine how many of the cells the signal has reached
sigma_D_all = [0.001 0.01 0.1];
nruns = 20;
t_max = 1000;
% number of cells that have turned at least one gene on (assuming background = (0,0))
num_cells_reached_all = cell(numel(sigma_D_all), nruns); 
folder = 'W:\staff-bulk\tnw\BN\HY\Shared\Yiteng\two_signals\moving_cells';
subfolder = 'subdomain_oscillations';
N = 15^2;

% Load data
for ii=1:numel(sigma_D_all)
    sigma_D = sigma_D_all(ii);
    for jj=1:nruns
        %fname_str = 'one_signal_sigma_D_0p001_t_out_1000-v1';
        fname_str = strrep(sprintf('two_signal_mult_sigma_D_%.3f_t_out_%d_period_Inf-v%d',...
            sigma_D, t_max, jj), '.' ,'p');
        fname = fullfile(folder, subfolder, strcat(fname_str, '.mat'));
        if exist(fname, 'file')==2
            disp(fname);
            load(fname)
            
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
%% Plot fractions over time
colors = {'b', 'r', 'g'};
h = figure;
p_temp = {};
hold on
for ii=1:numel(sigma_D_all)
    p_temp{ii} = plot(Inf, Inf, colors{ii});
    for jj=1:size(num_cells_reached_all, 2)
        frac_cells_reached = num_cells_reached_all{ii,jj}/N;
        plot(1:numel(frac_cells_reached), frac_cells_reached, colors{ii});
    end
end
legend([p_temp{1} p_temp{2} p_temp{3}], sprintfc('$\\sigma_D = %.3f$', sigma_D_all),...
    'Interpreter', 'latex');
%}