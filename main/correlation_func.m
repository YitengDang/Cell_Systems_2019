function [corr_func, dist_vals] = correlation_func(dist, cells)
% Computes the correlation function <Xi Xj> as a function of distance for a
% specific configuration
%% 
dist_vals = unique(round(dist, 5)); %unique distances; need to round to neglect numerical errors
corr_func = zeros(length(dist_vals), 1); % Correlation function
nterms = zeros(length(dist_vals), 1); % number of terms to average over for each value

%rows = {};
%cols = {};
for i=1:length(dist_vals)
    [row, col] = find(round(dist, 5) == dist_vals(i)); 
    %rows{end+1} = row;
    %cols{end+1} = col;
    nterms(i) = numel(row);
    %rows{i} = row; cols{i} = col;
    corr_func(i) = sum(cells(row).*cells(col))/numel(row) - mean(cells)^2; 
end

end