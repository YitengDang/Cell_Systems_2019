function[power_spec_dif] = calc_power_spec_dif(t_start,t_stop,cells_hist)
% For each time point of the simulation between t_start and t_stop two
% discrete spatial Fourier transformations are performed, 1 for each gene.
% This results in two power spectra. The absolute difference between these
% powerspectra gives an indication whether patterns can be present. In case
% of a wave the power spectra are identical and the difference between the
% power spectra would then be equal to 0. Structure that closely resamble
% waves would give values close to 0. However, homogeneous grids will also
% give a value of 0.

time_span = t_start+1:t_stop+1; % time 0 is in cells_hist{1}
power_spec_dif = zeros(1,length(time_span));

N = size(cells_hist{1},1); % gridsize
gz = sqrt(N);

%% Spatial Fourier transform to obtain power spectra
for t = time_span
    %cells = to_fourier_grid(round(cells_hist{t},12));
    cells = round(cells_hist{t},12);
    cells_norm = zeros(size(cells));
   
    std_1 = std(cells(:, 1));
    std_2 = std(cells(:, 2));
    
    % If grid is almost homogenous set to 0, otherwise oscillations in
    % aboslute difference due to position of the cells that are slightly
    % different
    if std_1 > 0.0001
        cells_norm(:, 1) = (cells(:, 1) - mean(cells(:, 1)))./std_1/gz;
    end
    
    if std_2 > 0.0001
        cells_norm(:, 2) = (cells(:, 2) - mean(cells(:, 2)))./std_2/gz;
    end
    
    cells_fft1 = fft2(reshape(cells_norm(:,1), gz, gz));
    cells_fft2 = fft2(reshape(cells_norm(:,2), gz, gz));
    
    power1 = round(abs(fftshift(cells_fft1)).^2/(N-1),12);
    power2 = round(abs(fftshift(cells_fft2)).^2/(N-1),12);
    
    % Absolute difference power spectra
    dif = calc_distance_to_pattern(power1,power2);
    
    power_spec_dif(t - t_start) = dif/2; % max dif could be 2 in principle, now 1.
end

end