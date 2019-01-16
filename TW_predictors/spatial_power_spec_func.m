function[power00,power10,power01,power11] = spatial_power_spec_func(cells)
    %%
    gz = sqrt(size(cells,1));
    N = size(cells,1);
    cell_states = zeros(size(cells,1),4);

    for i = 1:size(cells,1)
        state = translate_states(cells(i,:));
        cell_states(i,state) = 1;
    end
    %%
    %cell_states_norm = cell_states - mean(cell_states,1); OLD
    cell_states_norm = ((cell_states - mean(cell_states,1))./std(cell_states))/gz;

    cells_fft00 = fft2(reshape(cell_states_norm(:,1), gz, gz));
    %cells_fft1(1,1) = 0;
    cells_fft10 = fft2(reshape(cell_states_norm(:,2), gz, gz));
    %cells_fft2(1,1) = 0;
    cells_fft01 = fft2(reshape(cell_states_norm(:,3), gz, gz));
    cells_fft11 = fft2(reshape(cell_states_norm(:,4), gz, gz));

    power00 = abs(fftshift(cells_fft00)).^2/(N-1); % OLD: N
    power01 = abs(fftshift(cells_fft01)).^2/(N-1);
    power10 = abs(fftshift(cells_fft10)).^2/(N-1);
    power11= abs(fftshift(cells_fft11)).^2/(N-1);
end