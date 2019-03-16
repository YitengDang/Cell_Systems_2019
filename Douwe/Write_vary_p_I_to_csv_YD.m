%% Write data on initial p, I to csv file
clear all
close all
clc

%% Load data
nruns = 39;
folder = 'M:\tnw\bn\hy\Shared\Douwe\Recent Results\Vary_p_and_I';
fname_str = sprintf('analyzed_data_vary_p_I_nruns_%d', nruns);

load(fullfile(folder, fname_str), 'p1_all', 'p2_all', 'I1_ini_exact_all',...
    'I2_ini_exact_all', 'p_out_all', 'I_out_all', 'period_all', 't_out_all',...
    'trav_wave_all', 'trav_wave_2_all', 'names_all');

%% Store in vectors before printing to csv
p1_vec = [];
p2_vec = [];
I1_vec = [];
I2_vec = [];
TW_vec = [];
path_vec = {};
for p1_idx=1:numel(p1_all)
    for p2_idx=1:numel(p2_all)
        for I1_idx=1:size(I1_ini_exact_all, 3)
            for I2_idx=1:size(I1_ini_exact_all, 4)
                for run_idx=1:nruns
                    p1 = p1_all(p1_idx);
                    p2 = p2_all(p2_idx);
                    I1 = I1_ini_exact_all(p1_idx, p2_idx, I1_idx, I2_idx, run_idx);
                    I2 = I2_ini_exact_all(p1_idx, p2_idx, I1_idx, I2_idx, run_idx);
                    TW = trav_wave_2_all(p1_idx, p2_idx, I1_idx, I2_idx, run_idx);
                    
                    subfi1 = strrep(sprintf('p1_%.1f', p1), '.', 'p');
                    subfi2 = strrep(sprintf('p2_%.1f', p2), '.', 'p');
                    subfj1 = strrep(sprintf('I1_%.1f', I1), '.', 'p');
                    subfj2 = strrep(sprintf('I2_%.1f', I2), '.', 'p');
                    
                    fname_full = fullfile(folder, subfi1, subfi2, subfj1, subfj2,...
                        names_all{p1_idx, p2_idx, I1_idx, I2_idx, run_idx});
                    
                    p1_vec = [p1_vec, p1];
                    p2_vec = [p2_vec, p2];
                    I1_vec = [I1_vec, I1];
                    I2_vec = [I2_vec, I2];
                    TW_vec = [TW_vec, TW];
                    path_vec = [path_vec, fname_full];
                end
            end
        end
    end
end

%% Write to csv
save_folder = 'M:\tnw\bn\hy\Shared\Douwe\Recent Results\Vary_p_and_I';
M_all = {p1_vec',p2_vec',I1_vec',I2_vec',TW_vec',path_vec'};
var_name = {'p1', 'p2', 'I1', 'I2', 'TW', 'path'};

% Save File
% Produces a csv file with headers
cd(save_folder)
fid = fopen('Vary_I_p.csv','w');
fprintf(fid,'%s, %s, %s, %s, %s, %s\n',var_name{1,:})
for i = 1:length(p1_vec)
    fprintf(fid,'%f, %f, %f, %f, %f, %s\n',M_all{1}(i),M_all{2}(i),M_all{3}(i),M_all{4}(i),M_all{5}(i),M_all{6}{i})
end
fclose(fid)