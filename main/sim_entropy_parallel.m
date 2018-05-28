function [omega, omegak] = sim_entropy_parallel(dist, Son, K, a0, Rcell)
% Simulate the entropy by sampling several initial configurations by random

omega = 0; % initialize the number of eq. states to zeros
N = size(dist,1); % number of cells
n_smpl = 500;
omegak = zeros(N+1,1);

filename = fullfile(pwd,'data','nchoosek',sprintf('%d.mat', N));
if exist(filename, 'file') == 2
    tmp = load(filename);
    nck = tmp.nck;
else
    nck = zeros(N+1, 1);
    for n = 0:N
        nck(n+1) = nchoosek(N,n);
    end
    save(filename, 'nck');
end

for k = 0 : N
    if k == 0
        cells = zeros(N,1); % initialize all cells OFF
        [~, changed] = update_cells(cells, dist, Son, K, a0, Rcell);
        if ~changed
            omega = omega + 1;
            omegak(k+1) = 1;
        end
    else
        if k == N
            cells = ones(N,1); % initialize all cells ON
            [~, changed] = update_cells(cells, dist, Son, K, a0, Rcell);
            if ~changed
                omega = omega + 1;
                omegak(k+1) = 1;
            end
        else
            % Test n_smpl different sampling and see if it is eq.
            cnt_aux = zeros(n_smpl, 1);
            parfor i = 1:n_smpl
                rnd_idx = randperm(N,k);
                cells = zeros(N,1);
                cells(rnd_idx) = 1;
                [~, changed] = update_cells(cells, dist, Son, K, a0, Rcell);
                if ~changed
                    cnt_aux(i) = cnt_aux(i) +1;
                end
            end

            % the count of eq. states over the sampling number gives the
            % approximate probability of the number of eq. states with k ON
            % cells
            cnt = sum(cnt_aux);
            omegak(k+1) = cnt/n_smpl;
        end
    end
end
omegak = omegak.*nck;
omega = sum(sum(omegak));

end
