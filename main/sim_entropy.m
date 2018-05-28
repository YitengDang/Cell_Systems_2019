function [omega, omegak, nn] = sim_entropy(dist, Son, K, a0, Rcell)
% Simulate the entropy by sampling several initial configurations by random

omega = 0; % initialize the number of eq. states to zeros
N = size(dist,1); % number of cells
n_smpl = 500;
omegak = zeros(N+1,1);
nn = zeros(N+1,1);
% First neighbor distance
eps = 1e-5;
dist_vec = get_unique_distances(dist, eps);
dist1 = dist_vec(2);
idx = 1*(dist < dist1+eps & dist > dist1-eps);
for k = 0 : N
    if k == 0
        cells = zeros(N,1); % initialize all cells OFF
        [~, changed] = update_cells(cells, dist, Son, K, a0, Rcell);
        if ~changed
            omega = omega + 1;
            omegak(k+1) = 1;
            nn(k+1) = mean(idx*cells);
        end
    else
        if k == N
            cells = ones(N,1); % initialize all cells ON
            [~, changed] = update_cells(cells, dist, Son, K, a0, Rcell);
            if ~changed
                omega = omega + 1;
                omegak(k+1) = 1;
                nn(k+1) = mean(idx*cells);
            end
        else
            % Test n_smpl different sampling and see if it is eq.
            cnt = 0;
            for i = 1:n_smpl
                rnd_idx = randperm(N,k);
                cells = zeros(N,1);
                cells(rnd_idx) = 1;
                [~, changed] = update_cells(cells, dist, Son, K, a0, Rcell);
                if ~changed
                    cnt = cnt+1;
                    nn(k+1) = nn(k+1) + mean(idx*cells);
                end
            end

            % the count of eq. states over the sampling number gives the
            % approximate probability of the number of eq. states with k ON
            % cells
            if cnt > 0
                omegak(k+1) = cnt/n_smpl * nchoosek(N,k);
                nn(k+1) = nn(k+1)/cnt;
                omega = omega + omegak(k+1);
            end
        end
    end
end
omega = omega/(2^N);
end
