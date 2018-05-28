% Time evolution of a system with visualization of the dynamics without
% noise showing the count of nearest neighbors that are ON
% Few changes marked by 
%-------!!!!!!!!!!!!!!!!!!!!!--------------------------------------
close all
clear all
warning off
set(0, 'defaulttextinterpreter', 'latex');
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
K_loop = [2 5 10]; %2:15;
Con_loop = 5:5:40; %2:1:40;
p = [0.5 0.5];
Rcell = [0.2 0.2]*a0;

% Simulation parameters
nruns = 100; %100;
tmax = 1000;

% Communication matrix
Mcomm = [0 1; 1 0]; % type i reacts to type j iff M_ij=1
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

%-----------Calculate M-----------------
% cell types
cell_type = zeros(N,1);
idx1 = randperm(N,N1);
idx2 = setdiff(1:N, idx1);
cell_type(idx2) = 1;
% extra randomization of cell positions by requiring that cell_type also
% has a Moran's I value of around 0
[cell_type, ~] = generate_I_new(cell_type, 0, 0.01, dist, a0);
idx2 = find(cell_type);
idx1 = setdiff(1:N, idx2);

% Matrix of cell reading
M = zeros(size(dist));
    % same type
M(idx1, idx1) = Mcomm(1,1)*sinh(Rcell(1))./(a0*dist(idx1, idx1)).*exp(Rcell(1)-a0*dist(idx1, idx1));
M(sub2ind(size(M), idx1, idx1)) = Mcomm(1,1)*1; % 1: normal, 0: no self communication

M(idx2, idx2) = Mcomm(2,2)*sinh(Rcell(2))./(a0*dist(idx2, idx2)).*exp(Rcell(2)-a0*dist(idx2, idx2));
M(sub2ind(size(M), idx2, idx2)) = Mcomm(2,2)*1; % 1: normal, 0: no self communication
    % different type
M(idx1, idx2) = Mcomm(1,2)*Rcell(2)/(Rcell(1))*sinh(Rcell(1))./(a0*dist(idx1, idx2)).*exp(Rcell(2)-a0*dist(idx1, idx2)); % conc. cell of type 1 senses due to cells of type 2
M(idx2, idx1) = Mcomm(2,1)*Rcell(1)/(Rcell(2))*sinh(Rcell(2))./(a0*dist(idx2, idx1)).*exp(Rcell(1)-a0*dist(idx2, idx1));
    % self-communication
%M(sub2ind(size(M), 1:N, 1:N)) = 1; % 1: normal, 0: no self communication
%--------------------------------------------------------------------------
%% Load data
%subfolder = 'Mcomm_1_0_0_1_v2';
%path = fullfile('L:\BN\HY\Shared\Yiteng\Multicellularity\data\time_evolution_twocelltypes', subfolder);
path = 'L:\BN\HY\Shared\Yiteng\two_cell_types\data\time_evolution_twocelltypes\batch1\Mcomm_0_1_1_0';
listing = dir(path);
num_files = numel(listing)-2; %first two entries are not useful
straux = '(\d+)';
%K_fixed = [5 5];
%fpattern = strrep(...
%            sprintf('N1_%d_N2_%d_p1_%.2f_p2_%.2f_a0_%.2f_K1_%d_K2_%d_Con1_%s_Con2_%s_Rcell1_%.2f_Rcell2_%.2f_t%s-v%s', ...
%            N1, N2, p(1), p(2), a0, K_fixed(1), K_fixed(2), straux, straux, Rcell(1), Rcell(2), straux, straux), '.', 'p');      
fpattern = strrep(...
            sprintf('N1_%d_N2_%d_p1_%.2f_p2_%.2f_a0_%.2f_K1_%s_K2_%s_Con1_%s_Con2_%s_Rcell1_%.2f_Rcell2_%.2f_t%s-v%s', ...
            N1, N2, p(1), p(2), a0, straux, straux, straux, straux, Rcell(1), Rcell(2), straux, straux), '.', 'p');      
%fpattern = 'N1_113_N2_112_p1_0p50_p2_0p50_a0_1p50_K1_2_K2_2_Con1_(\d+)_Con2_(\d+)_Rcell1_0p30_Rcell2_0p30_t(\d+)-v(\d+)';        
count = 0;

I_final = zeros(numel(K_loop), numel(K_loop), numel(Con_loop), numel(Con_loop), nruns);
Non1_final = I_final;
Non2_final = I_final;
t_final = I_final;
theta_final = zeros(numel(K_loop), numel(K_loop), numel(Con_loop), numel(Con_loop), nruns, 3);
sim_count = zeros(numel(K_loop), numel(K_loop), numel(Con_loop), numel(Con_loop));
for i = 1:num_files
    filename = listing(i+2).name;
    % remove extension and do not include txt files
    [~,name,ext] = fileparts(filename);
    if strcmp(ext, '.mat')
        count = count + 1;
        if mod(count, 100)==0
            disp(count)
        end
        names{count} = name;
        [tokens, ~] = regexp(names{i}, fpattern,'tokens','match');
        if numel(tokens)==0
            disp('error');
            continue
        end
        % (1) tokens with K1, K2, Con1, Con2
        K1 = str2double(tokens{1}{1});
        K2 = str2double(tokens{1}{2});
        Con1 = str2double(tokens{1}{3});
        Con2 = str2double(tokens{1}{4});
        
        % (2) tokens with only Con1, Con2
        %K1 = K_fixed(1);
        %K2 = K_fixed(2);
        %Con1 = str2double(tokens{1}{1});
        %Con2 = str2double(tokens{1}{2});
        %-------!!!!!!!!!!!!!!!!!!!!!--------------------------------------
        %i1 = 1; i2 = 1;
        %-------!!!!!!!!!!!!!!!!!!!!!--------------------------------------
        i1 = find(K_loop==K1, 1);
        i2 = find(K_loop==K2, 1);
        i3 = find(Con_loop==Con1, 1);
        i4 = find(Con_loop==Con2, 1);
        sim_count(i1, i2, i3, i4) = sim_count(i1, i2, i3, i4) + 1;
        
        % store loaded data
        try 
            load(fullfile(path,filename), 'I', 'Non1', 'Non2', 't', ...
                'theta11', 'theta12', 'theta22', 'h_decrease');
        catch 
            fprintf('delete %s \n', filename)
            eval(sprintf('delete %s', fullfile(path,filename))); % delete corrupt files
        end
        % also delete files that have wrong data
        if (abs(I(end))>1)
            fprintf('delete %s \n', filename)
            eval(sprintf('delete %s', fullfile(path,filename))); % delete corrupt files
            continue; % NB skip 1 data point
        end
        
        count2 = sim_count(i1,i2,i3,i4);
        I_final(i1,i2,i3,i4,count2) = I(end);
        Non1_final(i1,i2,i3,i4,count2) = Non1(end);
        Non2_final(i1,i2,i3,i4,count2) = Non2(end);
        %theta_final(i1,i2,i3,i4,count2,:) = [theta11(end) theta12(end) theta22(end)];
        t_final(i1,i2,i3,i4,count2) = t;
    end
end
%% save statistics
qsave = 0;
if qsave 
    fname_str = strrep(sprintf('stats_out_N%d_a0_%.1f_N1_%d_nruns_%d_%s', N, a0, N1, nruns, subfolder), '.', 'p');
    path = 'L:\BN\HY\Shared\Yiteng\Multicellularity\data\time_evolution_twocelltypes';
    %path ='K:\bn\hy\Shared\Yiteng\Multicellularity\data_twocelltypes\time_evolution';
    fname = fullfile(path, fname_str);
    %fname = fullfile(pwd, 'data', 'time_evolution', fname_str);
    save(fname, 'Non1_final', 'Non2_final',...
        'theta_final', 'I_final', 'K_loop', 'Con_loop', 't_max');
end

%% Plot statistics
h1=figure(1);
histogram(squeeze(I_final(:)), -0.2:0.01:0.5);
xlabel('I');
%ylabel('Probability');
ylabel('Count');

h2=figure(2);
histogram(squeeze(Non1_final(:)), 0:N1);
xlabel('$$N_{ON, 1}$$');
ylabel('Count');
%ylabel('Probability');

h3=figure(3);
histogram(squeeze(Non2_final(:)), 0:N2);
xlabel('$$N_{ON, 2}$$');
%ylabel('Probability');
ylabel('Count');

h4=figure(4);
histogram(squeeze(t_final(:)), 0:tmax/20:tmax);
xlabel('$$t_{max}$$');
%ylabel('Probability');
ylabel('Count');
%% Calculate averages and plot
I_final_avg = zeros(numel(K_loop), numel(K_loop), numel(Con_loop), numel(Con_loop));
Non1_final_avg = I_final_avg;
Non2_final_avg = I_final_avg;
theta_final_avg = zeros(numel(K_loop), numel(K_loop), numel(Con_loop), numel(Con_loop), 3);
tmax_final_avg = I_final_avg;
for il1 = 1:numel(K_loop)
    for il2 = 1:numel(K_loop)
        for il3 = 1:numel(Con_loop)
            for il4 = 1:numel(Con_loop)
                I_final_avg(il1, il2, il3, il4) = mean(squeeze(I_final(il1,il2,il3,il4,:)));
                Non1_final_avg(il1, il2, il3, il4) = mean(squeeze(Non1_final(il1,il2,il3,il4,:)));
                Non2_final_avg(il1, il2, il3, il4) = mean(squeeze(Non2_final(il1,il2,il3,il4,:)));
                theta_final_avg(il1, il2, il3, il4, :) = mean(squeeze(theta_final(il1,il2,il3,il4,:,:)), 1);
                tmax_final_avg(il1, il2, il3, il4, :) = mean(squeeze(t_final(il1,il2,il3,il4,:)), 1);
            end
        end
    end
end
%}
%% Plot and save tmax
close all;
for k1=1:numel(K_loop)
    for k2=1:numel(K_loop)
        h=figure(numel(K_loop)*k1+k2); % figure index = 4, 5, ...
        imagesc(Con_loop, Con_loop, squeeze(tmax_final_avg(k1, k2, :, :))' );
        title(sprintf('$$K_1 = %d, K_2 = %d$$', K_loop(k1), K_loop(k2)));
        set(gca, 'YDir', 'normal', 'FontSize', 24);
        xlabel('$$C_{ON,1}$$')
        ylabel('$$C_{ON,2}$$')
        c=colorbar;
        ylabel(c, '$$\langle t_{final} \rangle$$', 'Interpreter', 'latex')
        caxis([0 tmax]);
        qsave = 1;
        if qsave 
            fname_str = strrep(sprintf('t_final_out_avg_N%d_a0_%.1f_N1_%d_K1_%d_K2_%d_nruns_%d_%s_v2',...
                N, a0, N1, K_loop(k1), K_loop(k2), nruns, subfolder), '.', 'p');
            %path = 'L:\BN\HY\Shared\Yiteng\Multicellularity\figures\time_evolution_twocelltypes';
            path = 'H:\My Documents\Multicellular automaton\two_cell_types\figures\time_evolution_twocelltypes\Mcomm_0_1_1_0';
            %path = 'L:\BN\HY\Shared\Yiteng\Multicellularity\data\time_evolution_twocelltypes';
            %path ='K:\bn\hy\Shared\Yiteng\Multicellularity\data_twocelltypes\time_evolution';
            fname = fullfile(path, fname_str);
            %fname = fullfile(pwd, 'data', 'time_evolution', fname_str);
            save_figure_pdf(h, 10, 8, fname);
        end
    end
end
%% Plot and save p1(tf)
close all;
for k1=1:numel(K_loop)
    for k2=1:numel(K_loop)
        h=figure(numel(K_loop)*k1+k2); % figure index = 4, 5, ...
        imagesc(Con_loop, Con_loop, squeeze(Non1_final_avg(k1, k2, :, :))'/N1 );
        title(sprintf('$$K_1 = %d, K_2 = %d$$', K_loop(k1), K_loop(k2)));
        set(gca, 'YDir', 'normal', 'FontSize', 24);
        xlabel('$$C_{ON,1}$$')
        ylabel('$$C_{ON,2}$$')
        c=colorbar;
        ylabel(c, '$$\langle N_{ON, 1}(t_f)/N_1 \rangle$$', 'Interpreter', 'latex')
        caxis([0 1]);
        qsave = 0;
        if qsave 
            fname_str = strrep(sprintf('Non1_out_avg_N%d_a0_%.1f_N1_%d_K1_%d_K2_%d_nruns_%d_%s',...
                N, a0, N1,  K_loop(k1), K_loop(k2), nruns, subfolder), '.', 'p');
            path = 'L:\BN\HY\Shared\Yiteng\Multicellularity\figures\time_evolution_twocelltypes';
            %path = 'L:\BN\HY\Shared\Yiteng\Multicellularity\data\time_evolution_twocelltypes';
            %path ='K:\bn\hy\Shared\Yiteng\Multicellularity\data_twocelltypes\time_evolution';
            fname = fullfile(path, fname_str);
            %fname = fullfile(pwd, 'data', 'time_evolution', fname_str);
            save_figure_pdf(h, 10, 8, fname);
        end
    end
end
%% Plot and save p2(tf)
close all;
for k1=1:numel(K_loop)
    for k2=1:numel(K_loop)
        h=figure(numel(K_loop)*k1+k2); % figure index = 4, 5, ...
        imagesc(Con_loop, Con_loop, squeeze(Non2_final_avg(k1, k2, :, :))'/N2 );
        title(sprintf('$$K_1 = %d, K_2 = %d$$', K_loop(k1), K_loop(k2)));
        set(gca, 'YDir', 'normal', 'FontSize', 24);
        xlabel('$$C_{ON,1}$$')
        ylabel('$$C_{ON,2}$$')
        c=colorbar;
        ylabel(c, '$$\langle N_{ON, 2}(t_f)/N_2 \rangle$$', 'Interpreter', 'latex')
        caxis([0 1]);
        qsave = 0;
        if qsave 
            fname_str = strrep(sprintf('Non2_out_avg_N%d_a0_%.1f_N1_%d_K1_%d_K2_%d_nruns_%d_%s',...
                N, a0, N1,  K_loop(k1), K_loop(k2), nruns, subfolder), '.', 'p');
            path = 'L:\BN\HY\Shared\Yiteng\Multicellularity\figures\time_evolution_twocelltypes';
            %path = 'L:\BN\HY\Shared\Yiteng\Multicellularity\data\time_evolution_twocelltypes';
            %path ='K:\bn\hy\Shared\Yiteng\Multicellularity\data_twocelltypes\time_evolution';
            fname = fullfile(path, fname_str);
            %fname = fullfile(pwd, 'data', 'time_evolution', fname_str);
            save_figure_pdf(h, 10, 8, fname);
        end
    end
end
%% Save theta^{11}
close all;
for k1=1:numel(K_loop)
    for k2=1:numel(K_loop)
        h=figure(numel(K_loop)*k1+k2); % figure index = 4, 5, ...
        imagesc(Con_loop, Con_loop, squeeze(theta_final_avg(k1, k2, :, :, 1))' );
        title(sprintf('$$K_1 = %d, K_2 = %d$$', K_loop(k1), K_loop(k2)));
        set(gca, 'YDir', 'normal', 'FontSize', 24);
        xlabel('$$C_{ON,1}$$')
        ylabel('$$C_{ON,2}$$')
        c=colorbar;
        %ylabel(c, '$$\langle N_{ON, 2}(t_f)/N_2 \rangle$$', 'Interpreter', 'latex')
        ylabel(c, '$$\langle \theta_{11}/(Nf^{11}) \rangle$$', 'Interpreter', 'latex')
        caxis([0 1]);
        qsave = 0;
        if qsave 
            fname_str = strrep(sprintf('Theta_11_out_avg_N%d_a0_%.1f_N1_%d_K1_%d_K2_%d_nruns_%d_%s',...
                N, a0, N1,  K_loop(k1), K_loop(k2), nruns, subfolder), '.', 'p');
            path = 'L:\BN\HY\Shared\Yiteng\Multicellularity\figures\time_evolution_twocelltypes';
            %path = 'L:\BN\HY\Shared\Yiteng\Multicellularity\data\time_evolution_twocelltypes';
            %path ='K:\bn\hy\Shared\Yiteng\Multicellularity\data_twocelltypes\time_evolution';
            fname = fullfile(path, fname_str);
            %fname = fullfile(pwd, 'data', 'time_evolution', fname_str);
            save_figure_pdf(h, 10, 8, fname);
        end
    end
end
%% Save theta^{12}
close all;
for k1=1:numel(K_loop)
    for k2=1:numel(K_loop)
        h=figure(numel(K_loop)*k1+k2); % figure index = 4, 5, ...
        imagesc(Con_loop, Con_loop, squeeze(theta_final_avg(k1, k2, :, :, 2))' );
        title(sprintf('$$K_1 = %d, K_2 = %d$$', K_loop(k1), K_loop(k2)));
        set(gca, 'YDir', 'normal', 'FontSize', 24);
        xlabel('$$C_{ON,1}$$')
        ylabel('$$C_{ON,2}$$')
        c=colorbar;
        %ylabel(c, '$$\langle N_{ON, 2}(t_f)/N_2 \rangle$$', 'Interpreter', 'latex')
        ylabel(c, '$$\langle \theta_{12}/(Nf^{12}) \rangle$$', 'Interpreter', 'latex')
        caxis([0 1]);
        qsave = 0;
        if qsave
            fname_str = strrep(sprintf('Theta_12_out_avg_N%d_a0_%.1f_N1_%d_K1_%d_K2_%d_nruns_%d_%s',...
                N, a0, N1, K_loop(k1), K_loop(k2), nruns, subfolder), '.', 'p');
            path = 'L:\BN\HY\Shared\Yiteng\Multicellularity\figures\time_evolution_twocelltypes';
            %path = 'L:\BN\HY\Shared\Yiteng\Multicellularity\data\time_evolution_twocelltypes';
            %path ='K:\bn\hy\Shared\Yiteng\Multicellularity\data_twocelltypes\time_evolution';
            fname = fullfile(path, fname_str);
            %fname = fullfile(pwd, 'data', 'time_evolution', fname_str);
            save_figure_pdf(h, 10, 8, fname);
        end
    end
end
%% Save theta^{22}
close all;
for k1=1:numel(K_loop)
    for k2=1:numel(K_loop)
        h=figure(numel(K_loop)*k1+k2); % figure index = 4, 5, ...
        imagesc(Con_loop, Con_loop, squeeze(theta_final_avg(k1, k2, :, :, 3))' );
        title(sprintf('$$K_1 = %d, K_2 = %d$$', K_loop(k1), K_loop(k2)));
        set(gca, 'YDir', 'normal', 'FontSize', 24);
        xlabel('$$C_{ON,1}$$')
        ylabel('$$C_{ON,2}$$')
        c=colorbar;
        %ylabel(c, '$$\langle N_{ON, 2}(t_f)/N_2 \rangle$$', 'Interpreter', 'latex')
        ylabel(c, '$$\langle \theta_{22}/(Nf^{22}) \rangle$$', 'Interpreter', 'latex')
        caxis([0 1]);
        qsave = 0;
        if qsave 
            fname_str = strrep(sprintf('Theta_22_out_avg_N%d_a0_%.1f_N1_%d_K1_%d_K2_%d_nruns_%d_%s',...
                N, a0, N1,  K_loop(k1), K_loop(k2), nruns, subfolder), '.', 'p');
            path = 'L:\BN\HY\Shared\Yiteng\Multicellularity\figures\time_evolution_twocelltypes';
            %path = 'L:\BN\HY\Shared\Yiteng\Multicellularity\data\time_evolution_twocelltypes';
            %path ='K:\bn\hy\Shared\Yiteng\Multicellularity\data_twocelltypes\time_evolution';
            fname = fullfile(path, fname_str);
            %fname = fullfile(pwd, 'data', 'time_evolution', fname_str);
            save_figure_pdf(h, 10, 8, fname);
        end
    end
end
%% Check 'wrong' data
%{
fpattern = strrep(...
            sprintf('N1_%d_N2_%d_p1_%.2f_p2_%.2f_a0_%.2f_K1_%d_K2_%d_Con1_%d_Con2_%d_Rcell1_%.2f_Rcell2_%.2f_t%s-v%s', ...
            N1, N2, p(1), p(2), a0, 5, 2, 10, 10, Rcell(1), Rcell(2), straux, straux), '.', 'p');      
I_temp = zeros(100,1);
index_temp = zeros(100,1);
count = 0;
for i=1:numel(files)    
	[tokens, ~] = regexp(files{i}, fpattern,'tokens','match');
    if numel(tokens)>0
        count = count+1;
        load(fullfile(path,files{i}), 'I', 'Non1', 'Non2', 't');
        I_temp(count) = I(end);
        index_temp(count) = i;
    end
end
%}