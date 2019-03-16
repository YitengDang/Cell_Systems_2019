%% Store in vectors before printing to csv
p1_vec = [];
p2_vec = [];
I1_vec = [];
I2_vec = [];
TW_vec = [];
path_vec = {};
%% Get all p1 folders (highest level) 
main_folder = 'W:\staff-groups\tnw\bn\hy\Shared\Douwe\Recent Results\Vary_p_and_I';
cd(main_folder)
p1_all = dir('p1_*');
p1_names = {p1_all.name};
%% Loop over all simulation files

for indx_p1 = 1:length(p1_names)
    current_p1 = sprintf('%s\\%s',main_folder,p1_names{indx_p1});
    
    % get current p2 folders
    cd(current_p1)
    p2_all = dir('p2_*');
    p2_names = {p2_all.name};
    
    for indx_p2 = 1:length(p2_names)
        current_p2 = sprintf('%s\\%s',current_p1,p2_names{indx_p2});
        
        % get current I1 folders
        cd(current_p2)
        I1_all = dir('I1_*');
        I1_names = {I1_all.name};
        
        for indx_I1 = 1:length(I1_names)
            current_I1 = sprintf('%s\\%s',current_p2,I1_names{indx_I1});
            
            % get current I2 folders
            cd(current_I1)
            I2_all = dir('I2_*');
            I2_names = {I2_all.name};
            
            for indx_I2 = 1:length(I2_names)
                current_I2 = sprintf('%s\\%s',current_I1,I2_names{indx_I2});
                
                % get all simulation files
                cd(current_I2)
                sim_all = dir('*.mat');
                sim_names = {sim_all.name};
                
                for indx_sim = 1:length(sim_names)
                    load(sim_names{indx_sim})
                    
                    % Test if TW
                    if t_onset == Inf
                        TW = 0;
                    else
                        TW = travelling_wave_test(cells_hist, save_consts_struct.a0, period, t_out,distances);
                        % This includes not entirely symmetric TWs with
                        % periods > gz
                    end
                    
                    % Get values
                    p1 = sum(cells_hist{1}(:,1))/length(cells_hist{1}(:,1));
                    p2 = sum(cells_hist{1}(:,2))/length(cells_hist{1}(:,2));
                    
                    Rcell = save_consts_struct.Rcell;
                    I1 = moranI(cells_hist{1}(:,1), distances,Rcell);
                    I2 = moranI(cells_hist{1}(:,2), distances,Rcell);
                    
                    full_path = sprintf('%s\\%s',current_I2,sim_names{indx_sim});
                    
                    % Store in vectors
                    p1_vec = [p1_vec, p1];
                    p2_vec = [p2_vec, p2];
                    I1_vec = [I1_vec, I1];
                    I2_vec = [I2_vec, I2];
                    TW_vec = [TW_vec, TW];
                    path_vec = [path_vec, full_path];
                    
                    
                end
            end
            
        end
        
    end
    
end
%% Write to csv
save_folder = 'D:\Documenten\Studie\Master\Internship Delft';

M_all = {p1_vec',p2_vec',I1_vec',I2_vec',TW_vec',path_vec'};
var_name = {'p1','p2','I1','I2','TW','path'};
%% Save File
% Produces a csv file with headers
cd(save_folder)
fid = fopen('Vary_I_p.csv','w');
fprintf(fid,'%s, %s, %s, %s, %s, %s\n',var_name{1,:})
for i = 1:length(p1_vec)
    fprintf(fid,'%f, %f, %f, %f, %f, %s\n',M_all{1}(i),M_all{2}(i),M_all{3}(i),M_all{4}(i),M_all{5}(i),M_all{6}{i})
end
fclose(fid)



