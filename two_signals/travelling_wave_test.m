function trav_wave = travelling_wave_test(cells_hist, a0, this_period, this_t_out)
    N = size(cells_hist{1}, 1);
    s = size(cells_hist{1}, 2);
    gz = sqrt(N);
    dist = init_dist_hex(gz, gz);
    
    trav_wave = 0;
    % (1) Period = gridsize
    if this_period == gz
        % (2) p(t), I(t) remain constant after onset periodicity 
        % calculate p, I
        p_last = zeros(this_period+1, 2);
        I_last = zeros(this_period+1, 2);

        for i1=1:this_period+1
            i2 = i1 + this_t_out-this_period;
            cells = cells_hist{i2};
            p_last(i1, :) = sum(cells, 1)/N;
            for j1=1:s
                [I_last(i1, j1), ~] = moranI(cells(:, j1), a0*dist);
            end
        end

        I_last = round(I_last, 5); % deal with rounding mistakes

        % check that it remains constant
        trav_wave = 1; % by default, assume it's a travelling wave
        p0 = p_last(1, :);
        I0 = I_last(1, :);
        for i1=1:this_period+1
            %disp( all(p_last(i1, :)==p0) )
            %disp( all(I_last(i1, :)==I0) )
            cond1 = all(p_last(i1, :)==p0);
            cond2 = all(I_last(i1, :)==I0);
            if ~cond1 || ~cond2
                trav_wave = 0; % not a travelling wave
                break
            end
        end
    end

end