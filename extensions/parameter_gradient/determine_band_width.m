function[bands_in_wave, number_of_waves] = determine_band_width(band_vec)
% Assume that if multiple waves, all waves equal number of bands.
band_nr = sum(band_vec);

% Determine how many bands are next to each other in a single wave
band_idx = find(band_vec == 1);
if ~isempty(band_idx)
    init_band = datasample(band_idx,1);
else
    bands_in_wave = 0;
    number_of_wave = 0;
    return
end

counter = 1;
% Check for neighbours on right (as positioned in band_vec)
stop = 0;
current_band = init_band; 

while stop ~= 1
    if current_band ~= length(band_vec)
        r_neighbour = current_band + 1;
    else
        r_neighbour = 1;
    end
    
    if band_vec(r_neighbour) == 1
        counter = counter + 1;
        current_band = r_neighbour;
    else
        stop = 1;
    end
end

% Check for neighbours on the left
stop = 0;
current_band = init_band; 

while stop ~= 1
    if current_band ~= 1
        l_neighbour = current_band - 1;
    else
        l_neighbour = length(band_vec);
    end
    
    if band_vec(l_neighbour) == 1
        counter = counter + 1;
        current_band = l_neighbour;
    else
        stop = 1;
    end
end

bands_in_wave = counter;

if mod((band_nr - bands_in_wave),bands_in_wave) == 0
    number_of_waves = band_nr / bands_in_wave;
else
    number_of_waves = NaN;
end



end

