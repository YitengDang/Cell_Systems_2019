function [wave_dist_vec, mean_dist_vec, background] = determine_wave_dist(t_start,t_stop,cells_hist)
% For each time point of the simulation between t_start and t_stop the wave
% distance is calculated. This is done in three steps. First, the four
% fourier spectra are calculated (1 for each state). Second, the background
% colour is determined, this is simply the most abundant state. Third, the
% pairwise distance between the three reamining spectra is summed. If there
% would be a wave this distance = 0. The more the grid resambles a wave
% pattern, the closer the value to 0.

%% States
% 00 = 1
% 10 = 2
% 01 = 3
% 11 = 4
time_span = t_start+1:t_stop+1; % time 0 is in cells_hist{1}
background = zeros(1,length(time_span));
wave_dist_vec = zeros(1,length(time_span));
%%

for t = time_span
    %cells_hist_four = to_fourier_grid(cells_hist{t});
    cells_hist_four = cells_hist{t};
    [power00,power10,power01,power11] = spatial_power_spec_func(cells_hist_four);
    power_states = {power00 power10 power01 power11};
    
    % calc majority state: background
    n_state_vec = [0 0 0 0];
    for s = 1:4
        for i = 1:size(cells_hist_four,1)
            state = translate_states(cells_hist_four(i,:));
            if state == s
                n_state_vec(s) = n_state_vec(s) + 1;
            end
        end
    end
        
    [~, maj_state] = max(n_state_vec);
    background(t-(time_span(1)-1)) = maj_state;
   
    wave_dist = 0;
    b = 1; % not both count 1-2 and 2-1 
    for s = 1:3
        if s ~= maj_state
            for s2 = b:4
                if s2 ~= maj_state
                    if s ~= s2
                        wave_dist = wave_dist + calc_distance_to_pattern(power_states{s},power_states{s2});
                    end
                end
            end
        end
        b = b + 1;
    end
    wave_dist_vec(t-(time_span(1)-1)) = wave_dist/6;
    % maximal difference between 2 power specs =2. 3 * 2 = 6.
end
mean_dist_vec = mean(wave_dist_vec);
std_dist_vec = std(wave_dist_vec);

%{
% Index holds 1 in almost wave if wave distance more than 3 std lower than
% the mean wave_distance.
almost_wave = zeros(1,length(wave_dist_vec));

for i = t_start:t_stop
    if wave_dist_vec(i) < mean_dist_vec - 3 * std_dist_vec
        almost_wave(i) = 1;
    end
end

% Get idecies (timepoints) where there almost_wave = 1
almost_wave = find(almost_wave == 1);
near_wave = {};

if length(unique(almost_wave)) > 1
    % If the indices in almost wave are consecetive timepoints, they are placed
    % in the same vector in near_wave. If they are isolated points, they will
    % be in a seperate vector in near wave.
    cur_wave = 1;
    while cur_wave < length(almost_wave)
        current_almost_wave = [almost_wave(cur_wave)];
        stop = 0;
        while stop ~= 1 && cur_wave < length(almost_wave)
            if almost_wave(cur_wave + 1) == almost_wave(cur_wave) + 1
                current_almost_wave = [current_almost_wave,almost_wave(cur_wave + 1)];
            else
                near_wave = [near_wave,current_almost_wave];
                stop = 1;
            end
            cur_wave = cur_wave + 1;
        end
    end
    
    if isempty(near_wave)
        near_wave = [near_wave,current_almost_wave];
    else
        if isequal(current_almost_wave,near_wave{length(near_wave)})
            near_wave = [near_wave,almost_wave(length(almost_wave))];
        else
            near_wave = [near_wave,current_almost_wave];
        end
    end
end

%}




end