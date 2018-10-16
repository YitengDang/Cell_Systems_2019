function [direction] = Determine_travel_direction(band_vec, band_vec2, orientation, bands_in_wave)
% inputs:
% band_vec -> vector with position of bands at time t
% band_vec2 -> vector with position of bands at time t + 1
if circshift(band_vec, 1*bands_in_wave) == band_vec2
    cshift = 1;
elseif circshift(band_vec, -1*bands_in_wave) == band_vec2
    cshift = -1;
else
    cshift = NaN;
end


if strcmp(orientation,'Horizontal') && cshift == 1
    direction = 'Up';
elseif strcmp(orientation,'Horizontal') && cshift == -1
    direction = 'Down';
elseif strcmp(orientation,'Vertical') && cshift == 1
    direction = 'Right';
elseif strcmp(orientation,'Vertical') && cshift == -1
    direction = 'Left';
else
    direction = NaN;
end

end