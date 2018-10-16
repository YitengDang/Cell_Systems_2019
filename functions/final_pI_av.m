%{
clear all
close all
% Load trajectory
path = 'L:\BN\HY\Shared\Yiteng\two_signals\sweep K12';
folder = '';
%name = 'two_signal_mult_N225_initiateI0_K12_24_t_out_6671_period_85-v1';
name = 'two_signal_mult_N225_initiateI0_K12_4_t_out_3_period_Inf-v11';
load(fullfile(path, folder, name), 'cells_hist');

this_period = Inf;
this_t_out = 6671;

% search for travelling waves
a0 = 1.5;
[p_av, pij_av, I_av] = final_pI_av_temp(cells_hist, a0, this_period, this_t_out);
%}
function [p_av, pij_av, I_av] = final_pI_av(cells_hist, a0, this_period, this_t_out)
    % Returns (average) final p, p_ij, I
    N = size(cells_hist{1}, 1);
    s = size(cells_hist{1}, 2);
    gz = sqrt(N);
    dist = init_dist_hex(gz, gz);
    
    if this_period==Inf
        % return final state values
        cells = cells_hist{end};
        p_av = sum(cells, 1)/N;
        I_av = zeros(s, 1);
        for j=1:s
            [I_av(j), ~] = moranI(cells(:, j), a0*dist);
        end
        
        % get p_ij
        cell_states = cells*[1; 2]; % conversion digital <-> binary
        % 0 <-> [0 0]
        % 1 <-> [1 0]
        % 2 <-> [0 1]
        % 3 <-> [1 1]
        pij_av = zeros(1, 4);
        for i=0:3
            pij_av(i+1) = sum(cell_states==i);
        end
        pij_av = pij_av/N;

    else
        % return averages
        p_last = zeros(this_period, 2);
        I_last = zeros(this_period, 2);
        p_ij_last = zeros(this_period, 4);
        
        for i1=1:this_period
            t1 = i1 + 1 + this_t_out - this_period; % last this_period time points
            cells = cells_hist{t1};
            p_last(i1, :) = sum(cells, 1)/N;
            for j1=1:s
                [I_last(i1, j1), ~] = moranI(cells(:, j1), a0*dist);
            end
            
            cell_states = cells*[1; 2];
            for i2=0:3
                p_ij_last(i1, i2+1) = sum(cell_states==i2);
            end
        end
        
        p_av = mean(p_last, 1);
        I_av = mean(I_last, 1);
        pij_av = mean(p_ij_last, 1)/N;
    end
end