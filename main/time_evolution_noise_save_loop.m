% Time evolution of a system with noise and without visualization
close all
clear all
warning off

% lattice parameters
gridsize = 11;
N = gridsize^2;
a0 = 0.5;
Rcell = 0.2*a0;

% simulation parameters
to_organize = 0;
n_run = 1;
tmax = 200;

% circuit parameters
%Con = 10;
%K = 10;
nrange = 0:N;
Klist = [10 10 14 19 16 18]; %(a0=1.5) [6 10 15 17 20]; % (a0=0.5) [10 10 14 19 16 18]; %[3];
Conlist = [5 21 16 14 8 6]; %(a0=1.5) [21 21 20 14 14]; % (a0=0.5) [5 21 16 14 8 6]; %[24]; 
%noise = 0.1*K; % noise calculated below

% initial condition
%p = 0.8;
%iniON = round(p*N);

% use hexagonal lattice
[dist, pos] = init_dist_hex(gridsize, gridsize);
dist_vec = a0*dist(1,:);
r = dist_vec(dist_vec>0); % exclude self influence
fN = sum(sinh(Rcell)*sum(exp(Rcell-r)./r)); % calculate signaling strength

% generate cell_type (0 case type 1, 1 case type 2)
cell_type = zeros(N,1); % all the same here

% calculate minimal noise and set noise
pminlist = 0:0.1:1; % list of p to examine minimal noise over
alphalist = zeros( numel(Klist), numel(pminlist) ); % list of minimum noise strengths
for i=1:numel(Klist)
    K = Klist(i);
    Con = Conlist(i);
    for j=1:numel(pminlist)
        pmin = pminlist(j);
        fprintf('pmin = %.1f \n', pmin);
        % pmin = 0.9;
        alpha_min = min(abs(Con+fN*(pmin*Con+1-pmin)-K), abs(1+fN*(pmin*Con+1-pmin)-K))/sqrt(N);
        alphalist(i,j) = 3*alpha_min;
        fprintf('alpha_min = %.2f \n', alphalist(i,j));
    end
end

alpha_list = max(alphalist, [], 2);
%}
%%
for i=1:numel(Klist)
    K = Klist(i);
    Con = Conlist(i);
    noise = alpha_list(i);
    for iniON = 0:N
        fprintf('n = %d \n', iniON);
        for tests = 1:n_run
            disp(tests)
            % initialize ON cells
            aux = round(to_organize*iniON);
            cells = zeros(N,1);
            cells(randperm(N,iniON-aux)) = 1;
            if aux > 0
                cells(find(cells==0,aux)) = 1;
            end

            % initialize figure
            %hin = figure(1);
            t = 0;
            I = [];
            Non = [];
            mom = [];
            cells_hist = {};
            %update_cell_figure(hin, pos, a0, cells, cell_type, t);
            cells_hist{end+1} = cells;
            Non(end+1) = sum(cells);
            I(end+1) = moranI(cells, a0*dist);
            [cells_out, changed, mom(end+1)] = update_cells_noise(cells, dist, Con, K, a0, Rcell, noise);
            for t = 1:tmax
                %k = waitforbuttonpress;
                %update_cell_figure(hin, pos, a0, cells_out, cell_type, t);
                cells_hist{end+1} = cells_out;
                I(end+1) = moranI(cells_out, a0*dist);
                Non(end+1) = sum(cells_out);
                cells = cells_out;
                [cells_out, changed, mom(end+1)] = update_cells_noise(cells, dist, Con, K, a0, Rcell, noise);
            end

            fname_str = strrep(sprintf('N%d_n%d_neq_%d_a0_%.1f_K_%d_Con_%d_noise_%.2f_t_%d', ...
                N, iniON, Non(end), a0, K, Con, noise, t), '.', 'p');
            i = 1;
            fname = fullfile(pwd, 'rebuttal', 'h_vector_field', 'data', ...
                strcat(fname_str,'-v',int2str(i),'.mat'));
            while exist(fname, 'file') == 2
                i=i+1;
                fname = fullfile(pwd, 'rebuttal', 'h_vector_field', 'data', ...
                    strcat(fname_str,'-v',int2str(i),'.mat'));
            end

            save(fname)
        end
    end
end