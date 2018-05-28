% Time evolution of a system without noise and without visualization

close all
clear all
warning off

%data_path = 'C:\Users\eduardopavinat\Dropbox\Matlab codes\data_onecelltype_entropy';

% lattice parameters
gridsize = 11;
N = gridsize^2;
a0 = 1.5;
Rcell = 0.2*a0;
% circuit parameters
Son = 21;
K = 6;
% initial conditions
p0 = 0.5;
I_target = 0.3;
iniON = round(p0*N);

% use hexagonal lattice
[dist, pos] = init_dist_hex(gridsize, gridsize);

dist_vec = a0*dist(1,:);
r = dist_vec(dist_vec>0); % exclude self influence
fN = sum(sinh(Rcell)*sum(exp(Rcell-r)./r)); % calculate signaling strength

% generate cell_type (0 case type 1, 1 case type 2)
cell_type = zeros(N,1); % all the same here

% --- initialize ON cells with given (p,I) ---
% search for saved cells, if not found, calculate new
path = fullfile(pwd, 'data', 'cells_ini_p_I');
straux = '(\d+)';
fpattern = strrep(sprintf('cells_a0_%.1f_p%.2f_I%.2f-v%s', a0, p0, I_target, straux), '.', 'p');

% Get all file names in the directory
listing = dir(path);
num_files = numel(listing)-2; %first two entries are not useful
count = 0;
for i = 1:num_files
    filename = listing(i+2).name;
    % remove extension and do not include txt files
    [~,name,ext] = fileparts(filename);
    if strcmp(ext, '.mat')
        count = count + 1;
        names{count} = name;
    end
end

% Load a random initial configuration with desired (p,I)
valid_files = zeros(num_files,1);
for i=1:num_files
    [tokens, ~] = regexp(names{i},fpattern,'tokens','match');
    if numel(tokens) > 0
        %disp(names{i}); % displays the file name
        valid_files(i) = 1;
    end    
end
valid_idx = find(valid_files==1);

%%
for tests = 1:sum(valid_files)
    % Start with random lattice with good $I$ each time
    disp(tests)
    this_idx = valid_idx(tests);
    disp(names{this_idx});
    
    load(fullfile(path,strcat(names{this_idx},'.mat')));  % load the data
    cells = cells_out;
    disp(strcat('I=', num2str(I_out)));
    
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
    [cells_out, changed, mom(end+1)] = update_cells(cells, dist, Son, K, a0, Rcell);
    while changed
        t = t+1;
        %k = waitforbuttonpress;
        %update_cell_figure(hin, pos, a0, cells_out, cell_type, t);
        cells_hist{end+1} = cells_out;
        I(end+1) = moranI(cells_out, a0*dist);
        Non(end+1) = sum(cells_out);
        cells = cells_out;
        [cells_out, changed, mom(end+1)] = update_cells(cells, dist, Son, K, a0, Rcell);
    end
    
    fname_str = strrep(sprintf('N%d_p0_%.2f_I0_%.2f_neq_%d_a0%.1f_K_%d_Son_%d_t_%d', ...
        N, iniON/N, I_out, Non(end), a0, K, Son, t), '.', 'p');
    i = 1;
    fname = fullfile(pwd, 'data','dynamics', 'no_noise_Itarget', ...
        strcat(fname_str,'-v',int2str(i),'.mat'));
    while exist(fname, 'file') == 2
        i=i+1;
        fname = fullfile(pwd, 'data', 'dynamics', 'no_noise_Itarget', ...
            strcat(fname_str,'-v',int2str(i),'.mat'));
    end
    
    save(fname)
end
