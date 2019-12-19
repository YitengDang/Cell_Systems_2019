% Time evolution of a system without noise and without visualization
close all
clear all
warning off

% Lattice parameters
gridsize = 11;
N = gridsize^2;
%a0 = 5;
%Rcell = 0.2*a0;
a0list = [4.5:0.1:4.9];
initialID = 'binaryrand'; % 'uniform';

% circuit parameters
K = 8;
Con = 16;
hill = 2; % Hill coefficient
prec = 8;
noise = 0;
%noiselist = 10.^[-1];

% simulation parameters
n_run = 100;
qsave = 1;

% use hexagonal lattice
[dist, pos] = init_dist_hex(gridsize, gridsize);

%%
for i=1:numel(a0list) %numel(noiselist)
    %noise = noiselist(i);
	a0 = a0list(i);
    Rcell = 0.2*a0;
    dist_vec = a0*dist(1,:);
    r = dist_vec(dist_vec>0); % exclude self influence
    fN = sum(sinh(Rcell)*sum(exp(Rcell-r)./r)); % calculate signaling strength

    % generate cell_type (0 case type 1, 1 case type 2)
    cell_type = zeros(N,1); % all the same here
    for tests = 1:n_run
        disp(tests)

        % initialize ON cells
        if strcmp(initialID, 'uniform')
            cells = rand(N, 1); %uniformly distributed
        elseif strcmp(initialID, 'fixedp')
            cells = p + sigma*randn(N, 1);
            cells(cells < 0) = 0;
            cells(cells > 1) = 1;
        elseif strcmp(initialID, 'fixedp_unif')
            cells = p*ones(N, 1);
        elseif strcmp(initialID, 'binary')
            cells = zeros(N, 1);
            cells(randperm(N, iniON)) = 1;
        elseif strcmp(initialID, 'binaryrand')
            cells = randi(2, N, 1) - 1;
        elseif strcmp(initialID, 'montecarlo')
            cells = init_random_cells_montecarlo(N, p);
        end

        t = 0;
        cells_hist = {};
        xmean = [];
        xstd = [];
        I = [];
        Theta = [];
        dY_mean = [];
        dY_std = [];
        
        cells_hist{end+1} = cells;
        xmean(end+1) = mean(cells); 
        xstd(end+1) = std(cells);
        [I(end+1), Theta(end+1)] = moranI(cells, a0*dist);
        [cells_out, changed, dY_mean(end+1), dY_std(end+1)] = ...
            update_cells_noise_hill(cells, dist, Con, K, a0, Rcell, noise, hill, prec);
        
        tmax = 1000;
        while t < tmax
            t = t+1;
            cells_hist{end+1} = cells_out;
            xmean(end+1) = mean(cells_out);
            xstd(end+1) = std(cells);
            [I(end+1), Theta(end+1)] = moranI(cells, a0*dist);
            cells = cells_out;
            [cells_out, changed, dY_mean(end+1), dY_std(end+1)] = ...
                update_cells_noise_hill(cells, dist, Con, K, a0, Rcell, noise, hill, prec);
        end

        % save file
        if qsave
            %fname_str = strrep(sprintf('N%d_a0_%.2f_K%d_Con%d_hill%.2f_t%d_xmeanf_%.2f_prec_%d_',...
            %    N, a0, K, Con, hill, t, xmean(end), prec), '.', 'p');
            %fname_str = strrep(sprintf('N%d_a0_%.2f_K%d_Con%d_hill%.2f_noise%.2f_prec_%d_tmax%d_%s',...
            %    N, a0, K, Con, hill, noise, prec, tmax, initialID), '.', 'p');
            fname_str = strrep(sprintf('N%d_a0_%.2f_Con%.2f_K%.2f_hill%.2f_%s',...
                N, a0, Con, K, hill, initialID), '.', 'p');
            i = 1;
            
            parent_folder = 'N:\tnw\BN\HY\Shared\Yiteng\one_signal_finite_Hill\dynamics';
            subfolder =  '2017-10-25_vs_a0_5p00_to_6p00_K8_Con16_hill2p00_binaryrand';
            save_dir = fullfile(parent_folder, subfolder); %'H:\My Documents\Multicellular automaton\temp';
            fname = fullfile(save_dir,...
                strcat(fname_str,'-v',int2str(i),'.mat'));
            while exist(fname, 'file') == 2
                i=i+1;
                fname = fullfile(save_dir,...
                    strcat(fname_str,'-v',int2str(i),'.mat'));
            end
            save(fname);
        end
    end
end