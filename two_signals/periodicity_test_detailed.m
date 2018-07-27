function [period, t_onset] = periodicity_test_detailed(cells_hist, t_check, period_ub)
    % for a confirmed periodic solution, checks within the last t_check+1
    % states when the first state revisit happened    
    
    % uses period_ub as an upper bound for the period obtained from periodicity_test_short
    t_out = numel(cells_hist)-1;
    for t1=t_out-t_check:t_out
        %disp(t1);
        cells_1 = cells_hist{t1+1};
        for t=t1-period_ub:t1-2 %check up to t1-2 (min. period = 2)
            %disp(t);
            cells_2 = cells_hist{t+1};
            if all(all(cells_2==cells_1))
                period = t1-t;
                t_onset = t;
                fprintf('t1=%d, period %d \n', t_onset, period);
                return
            end
        end
    end
    period = Inf;
    t_onset = Inf;
    warning('No period found in periodicity_test_detailed!')
    disp('Error: no period found');
end