function [xmean_out, xstd_out, t_out] = phenotype_sim2(dist, Con, K, a0, Rcell, hill)
    
nsim = 100; %number of simulations
N = size(dist,1); % number of cells
tmax = 1000; % maximum number of timesteps

% Initialize
xmean = zeros(1, nsim);
xstd = zeros(1, nsim);    
tvals = zeros(1, nsim);

for i=1:nsim
    % initialize ON cells
    %cells = rand(N, 1); %uniformly distributed
    cells = ones(N, 1); cells(1) = 0; %1 cell ON, rest OFF
    
    t = 0;   
    [cells_out, changed, ~] = update_cells_continuum(cells, dist, Con, K, a0, Rcell, hill);
    while changed && t<tmax
        t = t+1;
        cells = cells_out;
        [cells_out, changed, ~] = update_cells_continuum(cells, dist, Con, K, a0, Rcell, hill);
    end
    xmean(i) = mean(cells_out);
    xstd(i) = std(cells_out);
    tvals(i) = t;
end

%Mean values to export
xmean_out = mean(xmean);
xstd_out = mean(xstd);
t_out = mean(tvals);

end