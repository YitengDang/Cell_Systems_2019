function [trav_wave, trav_wave_2]  = travelling_wave_test(cells_hist, a0, this_period, this_t_out)
    % Predicts whether a given simulation is likely a travelling wave using
    % a set of criteria

    % Inputs:
    % (1) cells_hist: cell array with the simulation data
    % (2) a0: distance between cells
    % (3) this_period: period of the oscillation
    % (4) this_t_out: final time of the simulation 
    
    % Outputs:
    % In both cases, output = 1 if travelling wave is predicted, 0
    % otherwise
    % (1) trav_wave: strict criterion; p, I constant over one
    % period.
    % (2) trav_wave_2: loose criterion; p strictly constant.
    
    N = size(cells_hist{1}, 1);
    s = size(cells_hist{1}, 2);
    gz = sqrt(N);
    dist = init_dist_hex(gz, gz);
    
    trav_wave = 0;
    trav_wave_2 = 0;
    % Criterion 1: Period = (multiple of) gridsize
    if mod(this_period, gz)==0
        % Criterion 2: p(t), I(t) remain constant after onset periodicity 
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
        
        % deal with rounding mistakes in calculation of I
        I_last = round(I_last, 5);

        % check that the variables remain constant
        cond1 = size(unique(p_last, 'rows'), 1)==1;
        cond2 = size(unique(I_last, 'rows'), 1)==1;
        
        trav_wave = cond1 && cond2;
        trav_wave_2 = cond1;
        
        %{
        % old version
        %trav_wave = 1; % by default, assume it's a travelling wave
        %trav_wave_2 = 1;
        %p0 = p_last(1, :);
        %I0 = I_last(1, :);
        
        for i1=1:this_period+1
            %disp( all(p_last(i1, :)==p0) )
            %disp( all(I_last(i1, :)==I0) )
            cond1 = all(p_last(i1, :)==p0);
            cond2 = all(I_last(i1, :)==I0);
            if ~cond1
                trav_wave = 0; % not a travelling wave
                break
            end
        end
        %}
        % second test
        %{
        eps = 10^(-3);
        cond2b = all( std(I_last) < eps );
        if cond1 && cond2b
            trav_wave_2 = 1;
        end
        %}
    end

end