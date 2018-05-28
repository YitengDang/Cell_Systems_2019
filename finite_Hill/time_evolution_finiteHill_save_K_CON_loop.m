% Time evolution of a system without noise and without visualization
close all
clear all
warning off

% Lattice parameters
gridsize = 11;
N = gridsize^2;
%a0 = 5.5;
%Rcell = 0.2*a0;
initialID = 'uniform';
a0list = 3.5:0.5:6;

% circuit parameters
hill = 2; % Hill coefficient

% simulation parameters
n_run = 99;
qsave = 1;

% use hexagonal lattice
[dist, pos] = init_dist_hex(gridsize, gridsize);

%% Load K and Con values 
fname = fullfile(pwd, 'figures', 'finite_Hill_autonomy_KCon_map',...
    'Frac_autonomy_N121_a0_6p0_K5p00to15p00_Con10p00to25p00_hill2p00_uniform.mat');
load(fname, 'Klist_aut', 'Conlist_aut');
%%
for a0=a0list
    for j=1:numel(Klist_aut)
        K=Klist_aut(j);
        Con=Conlist_aut(j);
        fprintf('K=%d, Con=%d \n', K, Con);
        Rcell = 0.2*a0;
        dist_vec = a0*dist(1,:);
        r = dist_vec(dist_vec>0); % exclude self influence
        fN = sum(sinh(Rcell)*sum(exp(Rcell-r)./r)); % calculate signaling strength

        % generate cell_type (0 case type 1, 1 case type 2)
        cell_type = zeros(N,1); % all the same here
        %%
        for tests = 1:n_run
            %disp(tests)

            % initialize ON cells
            cells = rand(N, 1); %uniformly distributed
            % extend to other distributions

            t = 0;
            cells_hist = {};
            xmean = [];
            xstd = [];
            mom = [];
            I = [];

            cells_hist{end+1} = cells;
            xmean(end+1) = mean(cells); 
            xstd(end+1) = std(cells);
            I(end+1) = moranI(cells, a0*dist);
            [cells_out, changed, mom(end+1)] = update_cells_continuum(cells, dist, Con, K, a0, Rcell, hill);

            while changed && t < 10^4
                t = t+1;
                cells_hist{end+1} = cells_out;
                xmean(end+1) = mean(cells_out);
                xstd(end+1) = std(cells);
                I(end+1) = moranI(cells, a0*dist);
                cells = cells_out;
                [cells_out, changed, mom(end+1)] = update_cells_continuum(cells, dist, Con, K, a0, Rcell, hill);
            end

            % save file
            if qsave
                fname_str = strrep(sprintf('N%d_a0_%.2f_K%d_Con%d_hill%.2f_t%d_xmeanf_%.2f_',...
                    N, a0, K, Con, hill, t, xmean(end)), '.', 'p');
                i = 1;
                fname = fullfile(pwd, 'data','dynamics', 'finiteHill', '2017-09-08_2',...
                    strcat(fname_str,initialID,'-v',int2str(i),'.mat'));
                while exist(fname, 'file') == 2
                    i=i+1;
                    fname = fullfile(pwd, 'data', 'dynamics', 'finiteHill', '2017-09-08_2',...
                        strcat(fname_str,initialID,'-v',int2str(i),'.mat'));
                end
                avoidVariable = 'dist'; %avoid saving this variable
                save(fname, '-regexp', ['^(?!', avoidVariable,'$).'])
            end
        end
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