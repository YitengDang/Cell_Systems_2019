function [period, t_onset] = periodicity_test_short_reversed(cells_hist, decimals)
    % tests only if the last state has occured earlier by going back in
    % time (reversed direction)
        
    % decimals: optional rounding of cell states (for finite Hill)
    if nargin<2
        decimals = Inf;
    end
    
    t_out = numel(cells_hist) - 1;
    if decimals < Inf
        cells_current = round(cells_hist{end}, decimals);
    else
        cells_current = cells_hist{end};
    end
    
    for t=t_out-1:-1:0
        %disp(t);
        if decimals < Inf
            cells = round(cells_hist{t+1}, decimals);
        else
            cells = cells_hist{t+1};
        end
        if all(all(cells==cells_current))
            period = t_out-t;
            t_onset = t;
            fprintf('Periodicity test reversed: t_onset=%d, period: %d \n',...
                t_onset, period);
            return
        end
    end
    period = Inf;
    t_onset = Inf;
    %disp('no period found');
end