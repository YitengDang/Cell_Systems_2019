clear all
close all
set(0, 'defaulttextinterpreter', 'tex');

%% Load data
N = 225;
nruns = 39;
tmax = 10000;
a0 = 1.5;

p1_all = 0.1:0.2:0.9;
p2_all = 0.1:0.2:0.9;
%N1 = round(p1_all*N);
%N2 = round(p2_all*N);
I1_all = -0.1:0.1:0.7;
I2_all = -0.1:0.1:0.7;

%% Get all filenames
names_all = cell(numel(p1_all), numel(p2_all), numel(I1_all),...
 numel(I2_all), nruns);
counts_all = zeros(numel(p1_all), numel(p2_all), numel(I1_all),...
 numel(I2_all));
filepath = 'K:\bn\hy\Shared\Douwe\Recent Results\Vary_p_and_I';
for i1=1:numel(p1_all)
    p1 = p1_all(i1);
    subfi1 = strrep(sprintf('p1_%.1f', p1), '.', 'p');
    for i2=1:numel(p2_all)
        p2 = p2_all(i2);
        subfi2 = strrep(sprintf('p2_%.1f', p2), '.', 'p');
        for j1=1:numel(I1_all)
            I1 = I1_all(j1);
            subfj1 = strrep(sprintf('I1_%.1f', I1), '.', 'p');
            for j2=1:numel(I2_all)
                I2 = I2_all(j2);
                
                %{
                if I2==0
                    subfj2 = 'I2_0';
                else
                    subfj2 = strrep(sprintf('I2_%.1f', I2), '.', 'p');
                end
                %}
                subfj2 = strrep(sprintf('I2_%g', I2), '.', 'p');
                folder = fullfile(filepath, subfi1, subfi2, subfj1, subfj2);
                if exist(folder, 'dir')~=7
                    warning('folder does not exist, trying second naming convention');
                    subfj2 = strrep(sprintf('I2_%.1f', I2), '.', 'p');
                    folder = fullfile(filepath, subfi1, subfi2, subfj1, subfj2);
                    if exist(folder, 'dir')~=7 % try other folder
                        error('Folder does not exist!');
                    end
                end
                
                disp(folder);
                
                % Get all files in folder
                listing = dir(folder);
                num_files = numel(listing)-2; %first two entries are not useful
                count = 0;
                for i = 1:num_files
                    filename = listing(i+2).name;
                    % remove extension and do not include txt files
                    [~,name,ext] = fileparts(filename);
                    if strcmp(ext, '.mat')
                        count = count + 1;
                        names_all{i1, i2, j1, j2, count} = name;
                    end
                end
                counts_all(i1, i2, j1, j2) = count;
            end
        end
    end
end

%% Load data from files
filecount_all = zeros(numel(p1_all), numel(p2_all), numel(I1_all),...
 numel(I2_all));
I1_ini_exact_all = zeros(numel(p1_all), numel(p2_all), numel(I1_all),...
 numel(I2_all), nruns);
I2_ini_exact_all = I1_ini_exact_all;
trav_wave_all = I1_ini_exact_all;
trav_wave_2_all = I1_ini_exact_all;
period_all = I1_ini_exact_all;
t_out_all = I1_ini_exact_all;
p_out_all =  zeros(numel(p1_all), numel(p2_all), numel(I1_all),...
 numel(I2_all), nruns, 2);
I_out_all = p_out_all;
%%
for i1=1:numel(p1_all)
    p1 = p1_all(i1);
    subfi1 = strrep(sprintf('p1_%.1f', p1), '.', 'p');
    for i2=1:numel(p2_all)
        p2 = p2_all(i2);
        subfi2 = strrep(sprintf('p2_%.1f', p2), '.', 'p');
        for j1=1:numel(I1_all)
            I1 = I1_all(j1);
            subfj1 = strrep(sprintf('I1_%.1f', I1), '.', 'p');
            for j2=1:numel(I2_all)
                I2 = I2_all(j2);
                subfj2 = strrep(sprintf('I2_%.1f', I2), '.', 'p');                
                folder = fullfile(filepath, subfi1, subfi2, subfj1, subfj2);
                if exist(folder, 'dir')~=7
                    subfj2 = strrep(sprintf('I2_%g', I2), '.', 'p'); 
                    folder = fullfile(filepath, subfi1, subfi2, subfj1, subfj2);      
                end
                
                disp(folder);
                
                for k=1:nruns
                    fname_str = names_all{i1, i2, j1, j2, k};
                    fname = fullfile(folder, strcat(fname_str, '.mat'));
                    disp(fname);
                    
                    load(fname, 'cells_hist', 't_out', 'period', 'distances');
                    [trav_wave, trav_wave_2] = travelling_wave_test(cells_hist, a0,...
                        period, t_out, distances);
                    cells_ini = cells_hist{1};
                    
                    I1_ini_exact_all(i1, i2, j1, j2, k) = moranI(cells_ini(:,1),...
                        a0*distances);
                    I2_ini_exact_all(i1, i2, j1, j2, k) = moranI(cells_ini(:,2),...
                        a0*distances);
                    
                    trav_wave_all(i1, i2, j1, j2, k) = trav_wave;
                    trav_wave_2_all(i1, i2, j1, j2, k) = trav_wave_2;
                    period_all(i1, i2, j1, j2, k) = period;
                    t_out_all(i1, i2, j1, j2, k) = t_out;
                    
                    cells_out = cells_hist{end};
                    p_out_all(i1, i2, j1, j2, k, :) = mean(cells_out, 1);
                    I_out_all(i1, i2, j1, j2, k, 1) = moranI(cells_out(:,1), a0*distances);
                    I_out_all(i1, i2, j1, j2, k, 2) = moranI(cells_out(:,2), a0*distances);
                    
                    filecount_all(i1, i2, j1, j2) = filecount_all(i1, i2, j1, j2) + 1;
                end
            end
        end
    end
end

%% Save data
save_folder = 'K:\bn\hy\Shared\Douwe\Recent Results\Vary_p_and_I';
save_fname_str = sprintf('analyzed_data_vary_p_I_nruns_%d', nruns);
save( fullfile(save_folder, save_fname_str) );
%% Fix folder numbering issue
%
filepath = 'K:\bn\hy\Shared\Douwe\Recent Results\Vary_p_and_I';
for i1=1:numel(p1_all)
    p1 = p1_all(i1);
    %p1 = 0.3;
    subfi1 = strrep(sprintf('p1_%.1f', p1), '.', 'p');
    for i2=2:numel(p2_all)
        p2 = p2_all(i2);
        %p2 = 0.1;
        subfi2 = strrep(sprintf('p2_%.1f', p2), '.', 'p');
        for j1=1:numel(I1_all)
            %I1 = 0;
            I1 = I1_all(j1);
            subfj1 = strrep(sprintf('I1_%.1f', I1), '.', 'p');
            %subfj1 = 'I1_0p0';
            
            %
            subfj1_old = 'I1_0';
            folder_old = fullfile(filepath, subfi1, subfi2, subfj1_old);
            folder_new = fullfile(filepath, subfi1, subfi2, subfj1);
            %copyfile(folder_old, folder_new);
            if exist(folder_old, 'dir')==7 && exist(folder_new, 'dir')~=7
                disp(folder_old);
                copyfile(folder_old, folder_new);
            elseif exist(folder_old, 'dir')==7 && exist(folder_new, 'dir')==7
                fprintf('Removing directory %s \n', folder_old);
                rmdir(folder_old, 's');
            else 
                disp('Folder does not exist!')
            end
            %}
            
            %{
            for j2=1 %1:numel(I2_all)
                I2 = 0;
                subfj2_old = 'I2_0';
                subfj2 = 'I2_0p0';
                folder_old = fullfile(path, subfi1, subfi2, subfj1, subfj2_old);
                folder_new = fullfile(path, subfi1, subfi2, subfj1, subfj2);
                if exist(folder_old, 'dir')==7 && exist(folder_new, 'dir')~=7
                    disp(folder_old);
                    copyfile(folder_old, folder_new);
                elseif exist(folder_old, 'dir')==7 && exist(folder_new, 'dir')==7
                    fprintf('Removing directory %s \n', folder_old);
                    rmdir(folder_old, 's');
                else 
                    disp('Folder does not exist!')
                end
            end
            %}
        end
    end
end
%}