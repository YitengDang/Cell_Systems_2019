sim_folder = 'D:\Documenten\Studie\Master\Internship Delft\Robustness_patterns\Waves\Network_19\Initial\I_check\all_data\';

all_names = dir(sim_folder);
all_names = {all_names.name};
%%
I_all = -0.1:0.1:0.7;
a0 = 1.5;
gridsize = 15;

I1_vec = [];
I2_vec = [];
wave_vec = [];
t_wave_vec = [];
orientation_vec = {};
direction_vec = {};
nr_waves_vec = [];
nr_bands_vec = [];
band_colours_vec = {};
hamming_sync_vec = [];

for sim_nr = 1:100
    
    for I1 = 1:numel(I_all)
        
        for I2 = 1:numel(I_all)
            t_onset = NaN;
            
            % Load simulation
            I_ini_str = sprintf('_I_ini_%.2f_%.2f', I_all(I1), I_all(I2));
            file_name = strrep(sprintf('network_19_param_639_sim_%g_%s_',sim_nr,I_ini_str),'.','p');
            file_path = sprintf('%s%s%s',sim_folder,file_name,'*.mat');
            file_name = dir(file_path);
            file_name = file_name.name;
            load(fullfile(sim_folder, file_name))
            
            % Obtain wave properties in case of wave
            hamming_sync = calc_hamming_sync(cells_hist{1},distances,a0);
            if t_onset ~= Inf && t_out - t_onset >= gridsize
                wave = travelling_wave_test(cells_hist, a0, gridsize, t_onset + gridsize,distances);
            else
                wave = 0;
            end
            if wave
                [orientation,bands_in_wave,number_of_waves,diag_wrap,band_vec,wave_state,bended,wave_colours] = determine_wave_properties(cells_hist,t_onset);
                [~,~,~,~,band_vec2, ~,~,~] = determine_wave_properties(cells_hist,t_onset + 1,wave_state);
                if strcmp(orientation,'Horizontal') || strcmp(orientation,'Vertical')
                    direction = Determine_travel_direction(band_vec,band_vec2,orientation,bands_in_wave);
                elseif strcmp(orientation,'Diagonal')
                    direction_hor = Determine_travel_direction(band_vec{1},band_vec2{1},'Horizontal',sum(band_vec{1}));
                    direction_ver = Determine_travel_direction(band_vec{2},band_vec2{2},'Vertical',1);
                    direction = sprintf('%s & %s',direction_hor,direction_ver);
                else
                    direction = NaN;
                end
            else % in case no wave is detected
                orientation = NaN;
                bands_in_wave = NaN;
                number_of_waves = NaN;
                diag_wrap = NaN;
                band_vec = NaN;
                wave_colours = NaN;
                bended = NaN;
                band_vec2 = NaN;
                direction = NaN;
            end
            I1_vec = [I1_vec,I_all(I1)];
            I2_vec = [I2_vec,I_all(I2)];
            wave_vec = [wave_vec,wave];
            t_wave_vec = [t_wave_vec,t_onset];
            orientation_vec = [orientation_vec,orientation];
            direction_vec = [direction_vec,direction];
            nr_waves_vec = [nr_waves_vec,number_of_waves];
            nr_bands_vec = [nr_bands_vec,bands_in_wave];
            band_colours_vec = [band_colours_vec,wave_colours];
            hamming_sync_vec = [hamming_sync_vec,hamming_sync];
         
        end
    end
end


%% Heat matrix fraction wave
fraction_wave = sum(wave_vec)/length(wave_vec);
% Get fraction per unique I1 I2 combination
M_wave = [I2_vec',I1_vec',wave_vec',t_wave_vec']; % NOTE: gene 1 is 
% in col 2 and gene 2 in column 1, for ploting purposes
[C,ia,idx] = unique(M_wave(:,1:2),'rows');
M_wave_frac = [C,accumarray(idx,M_wave(:,3),[],@mean)];
M_t_wave = [C,accumarray(idx,M_wave(:,4),[],@mean)];
%%
I_val = size(unique(M_wave_frac(:,1)),1);
heat_matrix = zeros(I_val,I_val);
% Fill matrix
orig_index = 1;
for i = 1:I_val
    for j = 1:I_val
        heat_matrix(i,j) = M_wave_frac(orig_index,3);
        orig_index = orig_index + 1;
    end
end
%%
figure;
imagesc([M_t_wave(1,2) M_t_wave(size(M_t_wave,1),2)],[M_t_wave(1,1) M_t_wave(size(M_t_wave,1),1)],heat_matrix)
set(gca,'YDir','normal'); % y-axis form low to high
xlabel('I_0 gene 1')
ylabel('I_0 gene 2')
colorbar
caxis([0.0 1.0])
%% Boxplots 
I1_val = 0.8;
I2_val = 0.4;

I1_plot = M_wave(:,1) == I1_val;
I2_plot = M_wave(:,2) == I2_val; 
to_plot = I1_plot == 1 & I2_plot ==1;
M_plot = M_wave(to_plot,3);
%%
figure;
boxplot(M_plot)
hold on
scatter(ones(size(M_plot)).*(1+(rand(size(M_plot))-0.5)/10),M_plot,'k','filled')
set(gca,'xtick',[])
ylabel('Wave')
ylim([0.0 1.0])
title(sprintf('I1 = %g     I2 = %g',I1_val,I2_val));





%% Heatmap for t_wave
M_wave = [I2_vec(wave_vec == 1)',I1_vec(wave_vec == 1)',t_wave_vec(wave_vec == 1)'];
[C,ia,idx] = unique(M_wave(:,1:2),'rows');
M_t_wave = [C,accumarray(idx,M_wave(:,3),[],@mean)];
%%
I_val = size(unique(M_t_wave(:,1)),1);
heat_matrix = zeros(I_val,I_val);
% Fill matrix
orig_index = 1;
for i = 1:I_val
    for j = 1:I_val
        heat_matrix(i,j) = M_t_wave(orig_index,3);
        orig_index = orig_index + 1;
    end
end
%%
fig = figure;
imagesc([M_t_wave(1,2) M_t_wave(size(M_t_wave,1),2)],[M_t_wave(1,1) M_t_wave(size(M_t_wave,1),1)],heat_matrix)
set(gca,'YDir','normal'); % y-axis form low to high
xlabel('I_0 gene 1')
ylabel('I_0 gene 2')
colorbar
%% Hamming distance
M_ham= [round(hamming_sync_vec,1)',wave_vec']; 
[C,ia,idx] = unique(M_ham(:,1),'rows');
M_mean_ham = [C,accumarray(idx,M_ham(:,2),[],@mean)];
%%
figure
scatter(M_mean_ham(:,1),M_mean_ham(:,2))
%%
M_ham= [round(hamming_sync_vec(wave_vec == 1),1)',t_wave_vec(wave_vec == 1)']; 
[C,ia,idx] = unique(M_ham(:,1),'rows');
M_mean_ham = [C,accumarray(idx,M_ham(:,2),[],@mean)];

%% Save initial states sorted on wave and no wave
I_all = -0.1:0.1:0.9;
gridsize = sqrt(save_consts_struct.N);
a0 = save_consts_struct.a0;
Rcell = save_consts_struct.Rcell;
%%

for sim_nr = 1:8
    
    for I1 = 1:numel(I_all)
        
        for I2 = 1:numel(I_all)
            t_onset = NaN;
            
            % Load simulation
            I_ini_str = sprintf('_I_ini_%.2f_%.2f', I_all(I1), I_all(I2));
            file_name = strrep(sprintf('network_33_param_1890_sim_%g_%s_',sim_nr,I_ini_str),'.','p');
            file_path = sprintf('%s%s%s',sim_folder,file_name,'*.mat');
            file_name = dir(file_path);
            file_name = file_name.name;
            load(fullfile(sim_folder, file_name))
            
             if t_onset ~= Inf && t_out - t_onset >= gridsize
                wave = travelling_wave_test(cells_hist, a0, gridsize, t_onset + gridsize,distances);
            else
                wave = 0;
             end
            
             % Save image initial state
             %file_name = erase(file_name,".mat");
             if wave == 1
                 [orientation,~,~,~,~,~,~,~] = determine_wave_properties(cells_hist,t_onset);
                 if strcmp(orientation,'Horizontal')
                      save_folder = 'D:\Documenten\Studie\Master\Internship Delft\Robustness_patterns\Waves\Network_33\selected_network_33_param_1890_sim_2_t_out_1983_period_15\I_check\analysis\Init_states_wave\Horizontal\data';
                 elseif strcmp(orientation,'Vertical')
                      save_folder = 'D:\Documenten\Studie\Master\Internship Delft\Robustness_patterns\Waves\Network_33\selected_network_33_param_1890_sim_2_t_out_1983_period_15\I_check\analysis\Init_states_wave\Vertical\data';
                 elseif strcmp(orientation,'Diagonal')
                      save_folder = 'D:\Documenten\Studie\Master\Internship Delft\Robustness_patterns\Waves\Network_33\selected_network_33_param_1890_sim_2_t_out_1983_period_15\I_check\analysis\Init_states_wave\Diagonal\data';
                 else
                      save_folder = 'D:\Documenten\Studie\Master\Internship Delft\Robustness_patterns\Waves\Network_33\selected_network_33_param_1890_sim_2_t_out_1983_period_15\I_check\analysis\Init_states_wave\NaN';
                 end
                
             else
                 %continue
                 save_folder = 'D:\Documenten\Studie\Master\Internship Delft\Robustness_patterns\Waves\Network_33\selected_network_33_param_1890_sim_2_t_out_1983_period_15\I_check\analysis\Init_states_no_wave\data';
             end
            
             %{
             % plot fig
            h = figure;
            rcell = save_consts_struct.rcell;
            plot_handle = reset_cell_figure(h, positions, rcell);
            update_figure_periodic_scatter(plot_handle, cells_hist{1}, 0, 12, 1, a0, distances,Rcell);
            % save fig
            saveas(h, fullfile(save_folder, file_name), 'png')
             %}
             save(fullfile(save_folder, file_name),'cells_hist','distances','save_consts_struct','t_onset','t_out','positions','period')
        end
        %close all
    end
end
            




