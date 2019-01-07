function [period] = periodicity_test_short_forward(cells_hist, decimals)
    % tests only if the first state will occur again in the future by going
    % forward in time
        
    % decimals: optional rounding of cell states (for finite Hill)
    if nargin<2
        decimals = Inf;
    end
    
    t_out = numel(cells_hist) - 1;
    
    % get first state
    if decimals < Inf
        cells_current = round(cells_hist{1}, decimals);
    else
        cells_current = cells_hist{1};
    end
    
    % compare with forward states
    for t=1:t_out
        %disp(t);
        if decimals < Inf
            cells = round(cells_hist{t+1}, decimals);
        else
            cells = cells_hist{t+1};
        end
        if all(all(cells==cells_current))
            period = t;
            fprintf('Periodicity test forward: period: %d \n',...
            	period);
            return
        end
    end
    period = Inf;
    %disp('no period found');
end