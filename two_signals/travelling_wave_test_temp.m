clear all
close all
%% Load trajectory
path = 'L:\BN\HY\Shared\Yiteng\two_signals\sweep K12';
folder = '';
name = 'two_signal_mult_N225_initiateI0_K12_22_t_out_1090_period_15-v3';
load(fullfile(path, folder, name), 'cells_hist');

this_period = 15;
this_t_out = 1090;

% search for travelling waves
a0 = 1.5;
n_trav_wave = 0;
trav_wave = travelling_wave_test_temp(cells_hist, a0, this_period, this_t_out);

disp('Travelling wave?');
disp(trav_wave);
if trav_wave
    n_trav_wave = n_trav_wave + 1;
end

function trav_wave = travelling_wave_test_temp(cells_hist, a0, this_period, this_t_out)
    N = size(cells_hist{1}, 1);
    s = size(cells_hist{1}, 2);
    gz = sqrt(N);
    dist = init_dist_hex(gz, gz);
    
    trav_wave = 0;
    % (1) Period = gridsize
    if this_period == gz
        % (2) p(t), I(t) remain constant after onset periodicity 
        % calculate p, I
        p = zeros(numel(cells_hist), 2);
        I = zeros(numel(cells_hist), 1);

        for i1=1:numel(cells_hist)
            cells = cells_hist{i1};
            p(i1, :) = sum(cells, 1)/N;
            for j1=1:s
                [I(i1,j1), ~] = moranI(cells(:, j1), a0*dist);
            end
        end

        I = round(I, 5); % deal with rounding mistakes

        % check that it remains constant
        trav_wave = 1; % by default, assume it's a travelling wave
        p0 = p(this_t_out-this_period+1, :);
        I0 = I(this_t_out-this_period+1, :);
        for i1=this_t_out-this_period+2:this_t_out+1
            %disp( all(p(i1, :)==p0) )
            %disp( all(I(i1, :)==I0) )
            cond1 = all(p(i1, :)==p0);
            cond2 = all(I(i1, :)==I0);
            if ~cond1 || ~cond2
                trav_wave = 0; % not a travelling wave
                break
            end
        end
    end

end