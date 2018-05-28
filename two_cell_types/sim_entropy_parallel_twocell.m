function [omega, omegak] = sim_entropy_parallel_twocell(dist, Nv, Son, K, a0, Rcell)
% Simulate the entropy by sampling several initial configurations by random

N = size(dist,1); % number of cells
n_smpl = 500;
omegak = zeros(Nv+1);

nck_all = {[],[]};
i = 0;
for Nl = Nv
    i = i+1;
    filename = fullfile(pwd,'data','nchoosek',sprintf('%d.mat', Nl));
    if exist(filename, 'file') == 2
        tmp = load(filename);
        nck_all{i} = tmp.nck;
    else
        nck = zeros(Nl+1, 1);
        for n = 0:Nl
            nck(n+1) = nchoosek(Nl,n);
        end
        nck_all{i} = nck;
        save(filename, 'nck');
    end
end
[nck1, nck2] = meshgrid(nck_all{1}, nck_all{2});
N1 = Nv(1);
N2 = Nv(2);
for n1 = 0 : N1
    for n2 = 0 : N2
        % Test n_smpl different sampling and see if it is eq.
        cnt_aux = zeros(n_smpl, 1);
        K1 = K(1);
        K2 = K(2);
        Son1 = Son(1);
        Son2 = Son(2);
        parfor i = 1:n_smpl
            idx1 = randperm(N,N1);
            idx2 = setdiff(1:N, idx1);
            Kv = ones(N,1);
            Sonv = ones(N,1);
            Kv(idx1) = K1;
            Kv(idx2) = K2;
            Sonv(idx1) = Son1;
            Sonv(idx2) = Son2;
            cells = zeros(N,1);
            cells(idx1(randperm(N1,n1))) = 1;
            cells(idx2(randperm(N2,n2))) = 1;

            [~, changed] = update_cells(cells, dist, Sonv, Kv, a0, Rcell);
            if ~changed
                cnt_aux(i) = cnt_aux(i) +1;
            end
        end
        
        % the count of eq. states over the sampling number gives the
        % approximate probability of the number of eq. states with k ON
        % cells
        cnt = sum(cnt_aux);
        omegak(n1+1, n2+1) = cnt/n_smpl;
    end
end
omegak = omegak'.*nck1.*nck2;
omega = log(sum(sum(omegak)));

end
