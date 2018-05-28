% Time evolution of a system without noise and without visualization
% Start with binary cells
close all
clear all
warning off

% Lattice parameters
gridsize = 11;
N = gridsize^2;
a0 = 5;
Rcell = 0.2*a0;
initialID = 'binary';

% circuit parameters
K = 8;
Con = 16;
hill = 2; % Hill coefficient
prec = 8;

% simulation parameters
n_run = 1;
qsave = 1;

% use hexagonal lattice
[dist, pos] = init_dist_hex(gridsize, gridsize);
%%
iniONvals = 0:N;
for idx=1:numel(iniONvals)
    iniON = iniONvals(idx);
    dist_vec = a0*dist(1,:);
    r = dist_vec(dist_vec>0); % exclude self influence
    fN = sum(sinh(Rcell)*sum(exp(Rcell-r)./r)); % calculate signaling strength

    % generate cell_type (0 case type 1, 1 case type 2)
    cell_type = zeros(N,1); % all the same here
    for tests = 1:n_run
        disp(tests)

        % initialize ON cells
        cells = zeros(N,1);
        cells(randperm(N, iniON) ) = 1;

        t = 0;
        cells_hist = {};
        xmean = [];
        xstd = [];
        H = [];
        Theta = [];

        cells_hist{end+1} = cells;
        xmean(end+1) = mean(cells); 
        xstd(end+1) = std(cells);
        [~, Theta(end+1)] = moranI(cells, a0*dist);
        [cells_out, changed, H(end+1)] = update_cells_continuum(cells, dist, Con, K, a0, Rcell, hill, fN, prec);
        
        while changed && t < 10^4
            t = t+1;
            cells_hist{end+1} = cells_out;
            xmean(end+1) = mean(cells_out);
            xstd(end+1) = std(cells);
            [~, Theta(end+1)] = moranI(cells, a0*dist);
            cells = cells_out;
            [cells_out, changed, H(end+1)] = update_cells_continuum(cells, dist, Con, K, a0, Rcell, hill, fN, prec);
        end

        % save file
        if qsave
            fname_str = strrep(sprintf('N%d_n%d_a0_%.2f_K%d_Con%d_hill%.2f_prec_%d_',...
                N, iniON, a0, K, Con, hill, prec), '.', 'p');
            i = 1;
            fname = fullfile(pwd, 'data', '2017-10-31',...
                strcat(fname_str,initialID,'-v',int2str(i),'.mat'));
            while exist(fname, 'file') == 2
                i=i+1;
                fname = fullfile(pwd, 'data', '2017-10-31',...
                    strcat(fname_str,initialID,'-v',int2str(i),'.mat'));
            end
            save(fname);
        end
    end

    %{
    %% Plot <Xi>(t)
    figure();
    hold on
    plot(0:t, xmean, 'bo-');
    plot(0:t, xstd, 'ro-');

    %% Plot energy
    figure();
    plot(0:t, mom); %h = H/N
    %}
end