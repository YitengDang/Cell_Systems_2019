function [period, t_onset] = periodicity_test_short(cells_hist, save_consts_struct)
    % returns the first found period and the time of onset of the
    % periodicity
    N = save_consts_struct.N;
    gz = sqrt(N);
    t_out = numel(cells_hist)-1;
    % Scan over all initial frames
    %period = zeros(t_out+1, 1);
    for t1=0:t_out-2
        %disp(t1);
        %t1 = t_out-10; % initial frame
        cells_ini = cells_hist{t1+1};
        for t=t1+1:t_out
            %disp(t);
            cells = cells_hist{t+1};
            if all(all(cells==cells_ini))
                period = t-t1;
                t_onset = t1;
                fprintf('t1=%d, period %d \n', t1, t-t1);
                return
            end
        end
    end
    period = 0;
    disp('no period found');
end