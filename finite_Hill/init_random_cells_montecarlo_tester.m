% Confirms that mcsteps = 10^3 is a good choice (maximizes variance of
% cells)

clear all
N=100;
p=0.5;
ntrials = 100;
mcvals = 10.^(0:5);
var_vals = zeros(numel(mcvals), ntrials);
for i=1:numel(mcvals)
    for j=1:ntrials
        mcsteps = mcvals(i);
        cells = p*ones(N, 1);
        r = randi(N, mcsteps, 2);
        for k=1:mcsteps
            delta = (2*rand()-1)/2; % -0.5 <= delta <= 0.5
            cell1 = cells(r(k, 1)) + delta;
            cell2 = cells(r(k, 2)) - delta;
            if abs(2*cell1-1)<=1 && abs(2*cell2-1)<=1 % if both cells still in range [0, 1]
                cells(r(k, 1)) = cell1;
                cells(r(k, 2)) = cell2;
            end
        end
        var_vals(i, j) = var(cells);
    end
end

%%
varmean = mean(var_vals, 2);
varstd = zeros(numel(mcvals), 1);
for i=1:numel(mcvals)
    varstd(i) = std(var_vals(i,:))/sqrt(ntrials);
end

figure();
%plot(log10(mcvals), varmean);
errorbar(log10(mcvals), varmean, varstd);