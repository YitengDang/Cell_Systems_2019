clear all 
close all
%%
sim_folder = 'M:\tnw\bn\hy\Shared\Douwe\Recent Results\Initial_I\Network_19\all_data';
all_names = dir(sim_folder);
all_names = {all_names.name};
%%
I_all = -0.1:0.1:0.7;
I_dist_count = zeros( numel(I_all), numel(I_all), 100, 2 );  % get difference between target and actual I

a0 = 1.5;
gz = 15;
%rcell = 0.2;
%mcsteps = 0;
%[~, dist] = initial_cells_random_markov_periodic(gz, mcsteps, rcell);

wave_all = zeros( numel(I_all), numel(I_all), 100 );
wave_all2 = wave_all;
%{
I1_vec = [];
I2_vec = [];
wave_vec = [];
wave_vec2 = [];
t_wave_vec = [];
orientation_vec = {};
direction_vec = {};
nr_waves_vec = [];
nr_bands_vec = [];
band_colours_vec = {};
hamming_sync_vec = [];
%}
for sim_nr = 1:100
    disp(sim_nr);
    for I1 = 1:numel(I_all)
        
        for I2 = 1:numel(I_all)
            %t_onset = NaN;
            
            % Load simulation
            I_ini_str = sprintf('_I_ini_%.2f_%.2f', I_all(I1), I_all(I2));
            file_name = strrep(sprintf('network_19_param_639_sim_%g_%s_',sim_nr,I_ini_str),'.','p');
            file_path = sprintf('%s%s%s',sim_folder,file_name,'*.mat');
            file_path = sprintf('%s%s%s%s',sim_folder,'\', file_name,'*.mat');
            file_name_found = dir(file_path);
            
            %disp(file_name_found.name);
            load(fullfile(sim_folder, file_name_found.name))
            
            % check whether initial I matches target I (store distance)
            cells = cells_hist{1};
            I0 = [moranI(cells(:,1), a0*distances) moranI(cells(:,2), a0*distances)];
            %disp(I0);
            I0_target = [I_all(I1) I_all(I2)];
            I_dist_count(I1, I2, sim_nr, :) = I0 - I0_target;
            
            % Obtain wave properties in case of wave
            % hamming_sync = calc_hamming_sync(cells_hist{1},distances,a0);
            if mod(period, gz) == 0
                digits = 3;
                [wave, wave2] = travelling_wave_test(cells_hist, a0,...
                    period, t_out, distances, digits);
            else
                wave = 0;
                wave2 = 0;
            end
            fprintf('TW? wave = %d, wave2 = %d \n', wave, wave2);
            
            %{
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
            %}
            wave_all(I1, I2, sim_nr) = wave;
            wave_all2(I1, I2, sim_nr) = wave2;
            
            fprintf('Total? wave = %d, wave2 = %d \n',...
                sum(wave_all(I1, I2,:)), sum(wave_all2(I1, I2,:)));
            %{
            I1_vec = [I1_vec,I_all(I1)];
            I2_vec = [I2_vec,I_all(I2)];
            wave_vec = [wave_vec,wave];
            wave_vec2 = [wave_vec,wave2];
            %}
            % t_wave_vec = [t_wave_vec,t_onset];
            
            %{
            orientation_vec = [orientation_vec,orientation];
            direction_vec = [direction_vec,direction];
            nr_waves_vec = [nr_waves_vec,number_of_waves];
            nr_bands_vec = [nr_bands_vec,bands_in_wave];
            band_colours_vec = [band_colours_vec,wave_colours];
            hamming_sync_vec = [hamming_sync_vec,hamming_sync];
            %}
         
        end
    end
end

%% deviating I values
I_dist_temp = sum(abs(I_dist_count), 4);
I_dist_mean = mean(I_dist_temp, 3);

h=figure;
imagesc(I_all, I_all, I_dist_mean)
set(gca, 'YDir', 'normal')
colorbar;

%% make a selection
I_sel = 2:numel(I_all); % selection of I values to include
 % -> discard I=-0.1 values
I_sel_vals = I_all(I_sel);

cutoff = 0.05; % cut-off value
accepted_idx = I_dist_temp<cutoff;
n_accepted = sum(accepted_idx, 3);
frac_wave = zeros( numel(I_sel) ); 
frac_wave2 = frac_wave;

for ii1 = 1:numel(I_sel)
    for ii2 = 1:numel(I_sel)
        I1 = I_sel(ii1);
        I2 = I_sel(ii2);
        
        idx_temp = accepted_idx(I1, I2, :);
        n_accepted_temp = sum(idx_temp);

        n_wave = sum(wave_all(I1, I2, idx_temp), 3);
        frac_wave(ii1, ii2) = n_wave/n_accepted_temp;
        
        n_wave2 = sum(wave_all2(I1, I2, idx_temp), 3);
        frac_wave2(ii1, ii2) = n_wave2/n_accepted_temp;
        
    end
end

%% Save analyzed results
save_folder = 'H:\My Documents\Multicellular automaton\paper_2_draft\figures\originals_fig5';
fname_str = 'analyzed_trajectories_vs_I_network33_100runs';
save(fullfile(save_folder, fname_str), 'I_all', 'I_dist_count',...
    'wave_all', 'wave_all2', 'I_sel_vals', 'frac_wave', 'frac_wave2');

%% Heat map fraction wave
h=figure;
imagesc(I_all(I_sel), I_all(I_sel), frac_wave)
set(gca,'YDir','normal'); % y-axis form low to high
xlabel('Spatial order gene 1 (t=0)')
ylabel('Spatial order gene 2 (t=0)')
c=colorbar;
caxis([0.0 1.0])
ylabel(c, 'Fraction TW');
set(gca, 'XTick', I_all(I_sel));
set(gca, 'YTick', I_all(I_sel));
set(gca, 'FontSize', 24);
set(h, 'Units', 'Inches', 'Position', [0.1 0.1 10 8]);

% save figure
save_folder = 'H:\My Documents\Multicellular automaton\paper_2_draft\figures\originals_fig5';
qsave = 1;
fname_str = 'trajectories_vs_I_ini_frac_TW_strict_network19_100runs';
save_figure(h, 10, 8, fullfile(save_folder, fname_str),'.pdf', qsave);

%% Heat matrix fraction wave
%{
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
%}

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
gz = sqrt(save_consts_struct.N);
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
            
             if t_onset ~= Inf && t_out - t_onset >= gz
                wave = travelling_wave_test(cells_hist, a0, gz, t_onset + gz,distances);
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
            




