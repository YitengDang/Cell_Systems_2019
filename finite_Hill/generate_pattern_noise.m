function [cells_out] = generate_pattern_noise(cells_in, noise)
    %clear all
    %cells_in = zeros(100, 1);
    %cells_in(randperm(100, 50))=1;
    %noise = 1;
    N = numel(cells_in);
    cells_out = cells_in + noise*randn(N, 1);
    idx = find((cells_out > 1) + (cells_out < 0));
    while ~isempty(idx)
        cells_out(idx) = cells_in(idx) + noise*randn(numel(idx), 1);
        idx = find((cells_out > 1) + (cells_out < 0));
    end
    %    
    %figure();
    %plot(cells_in, cells_out, 'x');
end