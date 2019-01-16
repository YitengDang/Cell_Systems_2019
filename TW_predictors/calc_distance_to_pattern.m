function[distance] = calc_distance_to_pattern(current_spec,pattern_spec)
difference = abs(current_spec - pattern_spec);
distance = sum(sum(difference));
end