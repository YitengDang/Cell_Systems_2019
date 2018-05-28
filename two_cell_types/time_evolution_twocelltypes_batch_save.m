% Batch simulations of a large collection of trajectories
close all
clear all
warning off

%--------------------------------------------------------------------------
% Set parameters of the system
gridsize = 15;
N = gridsize^2;
frac_type1 = 0.5;
N1 = round(frac_type1*N);
N2 = N - N1;
a0 = 1.5;

% Parameters 
%Con = [20 16]; %column j = type j 
%K = [6 8];
p = [0.5 0.5];
Rcell = [0.2 0.2]*a0;

% Simulation parameters
nruns = 100; %100;
tmax = 1000;

% Communication matrix
Mcomm = [1 1; 1 1]; % type i reacts to type j iff M_ij=1
Mcomm = logical(Mcomm);
subfolder = sprintf('Mcomm_%d_%d_%d_%d', Mcomm(1,1), Mcomm(1,2), Mcomm(2,1), Mcomm(2,2));

% Check input
if ~all(size(Mcomm)==[2 2])
    disp('Wrong communication matrix size!');
    return
end
for idx=1:numel(Mcomm)
    if ~(Mcomm(idx)==0 || Mcomm(idx)==1)
        disp('Wrong communication matrix entries!');
        return
    end
end

% Distance and position
[dist, ~] = init_dist_hex(gridsize, gridsize);

% cell types
%cell_type = zeros(N,1);
%idx1 = randperm(N,N1);
%idx2 = setdiff(1:N, idx1);
%cell_type(idx2) = 1;
% extra randomization of cell positions by requiring that cell_type also
% has a Moran's I value of around 0
%[cell_type, ~] = generate_I_new(cell_type, 0, 0.01, dist, a0);
%idx2 = find(cell_type);
%idx1 = setdiff(1:N, idx2);
%--------------------------------------------------------------------------
%% Main loop
K_loop = 10; %[2 5 10]; %2:15;
Con_loop = 5:5:40; %2:1:40;
dt = 1.4; % estimated time per li5
t_est = nruns*numel(K_loop)^2*numel(Con_loop)^2*dt;
fprintf('Estimated running time: %.f s = %.2f min = %.2f h \n', t_est, t_est/60, t_est/3600);
I_final = zeros(numel(K_loop), numel(K_loop), numel(Con_loop), numel(Con_loop), nruns);
Non1_final = I_final;
Non2_final = I_final;
t_max = I_final;
h_decrease = I_final;
theta_final = zeros(numel(K_loop), numel(K_loop), numel(Con_loop), numel(Con_loop), nruns, 3);
%K = [5 5];
%Con = [35 5]; %zeros(2,1);
savepath = fullfile('L:\HY\Shared\Yiteng\Multicellularity\data\time_evolution_twocelltypes', subfolder);
%path = fullfile('L:\BN\HY\Shared\Yiteng\Multicellularity\data\time_evolution_twocelltypes\temp');
%%
for li1 = 1:numel(K_loop)
    K(1) = K_loop(li1);
    for li2 = 1:numel(K_loop)
        K(2) = K_loop(li2);
        for li3 = 1:numel(Con_loop)
            Con(1) = Con_loop(li3);
            for li4 = 1:numel(Con_loop)
                Con(2) = Con_loop(li4);
                fprintf('K1 = %d, K2 = %d, Con1 = %d, Con2 = %d \n', K(1), K(2), Con(1), Con(2));
                %fprintf('li = %d, %d, %d, %d', li1, li2, li3, li4);
                for li5=1:nruns
                    disp(li5);
                    [Non1_final(li1,li2,li3,li4,li5), Non2_final(li1,li2,li3,li4,li5),...
                        I_final(li1,li2,li3,li4,li5), t_max(li1,li2,li3,li4,li5),...
                        h_decrease(li1,li2,li3,li4,li5), theta_final(li1,li2,li3,li4,li5,:)] = ...
                        time_evolution_twocelltypes_full(N, dist,...
                        a0, Con, K, Rcell, frac_type1, Mcomm, p, tmax, savepath);
                    %[~] = time_evolution_twocelltypes_full(N, dist,...
                    %   a0, Con, K, Rcell, frac_type1, p, M);
                end
                %}
            end
        end
    end
end
%}
%% Get missing data points
%
sim_count = zeros(numel(K_loop), numel(K_loop), numel(Con_loop), numel(Con_loop));
%path = 'L:\HY\Shared\Yiteng\Multicellularity\data\time_evolution_twocelltypes\Mcomm_0_1_1_0';
listing = dir(savepath);
num_files = numel(listing)-2; %first two entries are not useful
count = 0;
straux = '(\d+)';
fpattern = strrep(...
            sprintf('N1_%d_N2_%d_p1_%.2f_p2_%.2f_a0_%.2f_K1_%s_K2_%s_Con1_%s_Con2_%s_Rcell1_%.2f_Rcell2_%.2f_t%s-v%s', ...
            N1, N2, p(1), p(2), a0, straux, straux, straux, straux, Rcell(1), Rcell(2), straux, straux), '.', 'p');      
for i = 1:num_files
    filename = listing(i+2).name;
    % remove extension and do not include txt files
    [~,name,ext] = fileparts(filename);
    if strcmp(ext, '.mat')
        count = count + 1;
        names{count} = name;
        [tokens, ~] = regexp(names{i}, fpattern,'tokens','match');
        if numel(tokens)==0
            disp('error');
            continue
        end
        K1 = str2double(tokens{1}{1});
        K2 = str2double(tokens{1}{2});
        Con1 = str2double(tokens{1}{3});
        Con2 = str2double(tokens{1}{4});
        i1 = find(K_loop==K1, 1);
        i2 = find(K_loop==K2, 1);
        i3 = find(Con_loop==Con1, 1);
        i4 = find(Con_loop==Con2, 1);
        sim_count(i1,i2,i3,i4) = sim_count(i1,i2,i3,i4)+1;
    end
end

% Redo missing points
for li1 = 1:numel(K_loop)
    K(1) = K_loop(li1);
    for li2 = 1:numel(K_loop)
        K(2) = K_loop(li2);
        for li3 = 1:numel(Con_loop)
            Con(1) = Con_loop(li3);
            for li4 = 1:numel(Con_loop)
                Con(2) = Con_loop(li4);
                fprintf('K1 = %d, K2 = %d, Con1 = %d, Con2 = %d \n', K(1), K(2), Con(1), Con(2));
                runstodo = nruns-sim_count(li1,li2,li3,li4);
                fprintf('Runs to do = %d \n', runstodo);
                for li5=1:runstodo
                    disp(li5);
                    %[Non1_final, Non2_final, I_final, tmax_out, h_decrease, theta_final] = time_evolution_twocelltypes_full(N, dist,...
                    %    a0, Con, K, Rcell, frac_type1, Mcomm, p, tmax, savepath);
                    [Non1_final(li1,li2,li3,li4,li5), Non2_final(li1,li2,li3,li4,li5),...
                        I_final(li1,li2,li3,li4,li5), t_max(li1,li2,li3,li4,li5),...
                        h_decrease(li1,li2,li3,li4,li5), theta_final(li1,li2,li3,li4,li5,:)] ...
                        = time_evolution_twocelltypes_full(N, dist,...
                       a0, Con, K, Rcell, frac_type1, Mcomm, p, tmax, savepath);
                    %[~] = time_evolution_twocelltypes_full(N, dist,...
                    %   a0, Con, K, Rcell, frac_type1, p, M);
                end
            end
        end
    end
end
%}

% save statistics
qsave = 1;
if qsave
    fname_str = strrep(sprintf('stats_out_N%d_a0_%.1f_N1_%d_nruns_%d_%s', N, a0, N1, nruns, subfolder), '.', 'p');
    savepath = 'L:\HY\Shared\Yiteng\Multicellularity\data\time_evolution_twocelltypes';
    %path ='K:\bn\hy\Shared\Yiteng\Multicellularity\data_twocelltypes\time_evolution';
    fname = fullfile(savepath, fname_str);
    %fname = fullfile(pwd, 'data', 'time_evolution', fname_str);
    save(fname, 'Non1_final', 'Non2_final', 'I_final', 'h_decrease', ...
        'theta_final', 'K_loop', 'Con_loop');
end

%% run with random p
%{
p_final = zeros(numel(K_loop), numel(K_loop), numel(Con_loop), numel(Con_loop), nruns);
I_final = p_final;
K = zeros(2,1);
Con = zeros(2,1);
for li1 = 1:numel(K_loop)
    K(1) = K_loop(li1);
    for li2 = 1:numel(K_loop)
        K(2) = K_loop(li2);
        for li3 = 1:numel(Con_loop)
            Con(1) = Con_loop(li3);
            for li4 = 1:numel(Con_loop)
                Con(2) = Con_loop(li4);
                fprintf('K1 = %d, K2 = %d, Con1 = %d, Con2 = %d \n', K(1), K(2), Con(1), Con(2));
                for li5=1:nruns
                    disp(li5);
                    p = rand(1, 2);
                    [Non1_final(li1,li2,li3,li4,li5), Non2_final(li1,li2,li3,li4,li5),...
                        I_final(li1,li2,li3,li4,li5)] = time_evolution_twocelltypes_full(N, dist,...
                       a0, Con, K, Rcell, frac_type1, p, M);
                    %[~] = time_evolution_twocelltypes_full(N, dist,...
                    %   a0, Con, K, Rcell, frac_type1, p, M);
                end
            end
        end
    end
end

% save statistics
qsave = 1;
if qsave 
    fname_str = strrep(sprintf('stats_out_randp_N%d_a0_%.1f_N1_%d_nruns_%d', N, a0, N1, nruns), '.', 'p');
    path = 'L:\HY\Shared\Yiteng\Multicellularity\data\time_evolution_twocelltypes\randp';
    %path ='K:\bn\hy\Shared\Yiteng\Multicellularity\data_twocelltypes\time_evolution';
    fname = fullfile(path, fname_str);
    %fname = fullfile(pwd, 'data', 'time_evolution', fname_str);
    save(fname, 'Non1_final', 'Non2_final', 'I_final');
end
%}