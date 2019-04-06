load_folder = 'D:\Documenten\Studie\Master\Internship Delft\Robustness_patterns\Waves\Network_33\selected_network_33_param_1890_sim_2_t_out_1983_period_15\Hill';

cd(load_folder)
file_names = dir('*.mat');
file_names = {file_names.name};
%%
per = 15;
wave_vec_simple = zeros(1,length(file_names));
wave_vec_fourier = zeros(1,length(file_names));
hill_vec = zeros(1,length(file_names));

t_out_vec = [];
t_onset_vec = [];
orientation_vec = {};
direction_vec = {};
nr_waves_vec = [];
nr_bands_vec = [];
bended_vec = [];
diag_wrap_vec = {};
wave_colours_vec = {};
I1_vec = {};
I2_vec = {};

for k = 1:length(file_names)
    load(file_names{k})
    
    hill_vec(k) = save_consts_struct.hill;
    % Only check if wave in case of periodicity
    wave_simple = 0;
    if t_onset ~= Inf
        
        for t = 1:length(cells_hist)
            cells_hist{t} = round(cells_hist{t},12);
        end
      
        wave_simple = travelling_wave_test(cells_hist, save_consts_struct.a0, per, t_out,distances);
      
        %wave_four = wave_test_fourier(cells_hist,per,t_out);
        
        wave_vec_simple(k) = wave_simple;
        %wave_vec_fourier(k) = wave_four;
    end
    %{
    if wave_simple == 1
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
    else % no wave detected
        orientation = -1;
        bands_in_wave = -1;
        number_of_waves = -1;
        diag_wrap = -1;
        band_vec = -1;
        wave_colours = -1;
        bended = -1;
        band_vec2 = -1;
        direction = -1;
        
    end
    %}
    t_out_vec = [t_out_vec,t_out];
    t_onset_vec = [t_onset_vec,t_onset];
    
    I1_values = zeros(0,length(cells_hist));
    I2_values = zeros(0,length(cells_hist));
    for i = 1:length(cells_hist)
        I1_values(i) = moranI(cells_hist{i}(:,1), distances,save_consts_struct.Rcell);
        I2_values(i) = moranI(cells_hist{i}(:,2), distances,save_consts_struct.Rcell);
    end
    I1_vec = [I1_vec,I1_values];
    I2_vec = [I2_vec,I2_values];
    %orientation_vec = [orientation_vec,orientation];
    %direction_vec = [direction_vec, direction];
    %nr_waves_vec = [nr_waves_vec,number_of_waves];
    %nr_bands_vec = [nr_bands_vec,bands_in_wave];
    %bended_vec = [bended_vec,bended];
    %diag_wrap_vec = [diag_wrap_vec,diag_wrap];
    %wave_colours_vec = [wave_colours_vec,wave_colours];
    
    fprintf(sprintf('k %g from %g \n',k,length(file_names)));
end
%%
save('Correct_Overview_rounded_p_2_dec.mat')
%% Check for what Hill coefficients NaN appears for the genes
grid_nan = zeros(1,5); % for every hill-coef how many cell grids contain one or more nan values [5 10 50 100 inf]
for k = 1:length(file_names)
    load(file_names{k})
    fprintf(sprintf('k %g from %g \n',k,length(file_names)));
    
    for t = 1:length(cells_hist)
       grid = cells_hist{t};
       nan_sum = sum(sum(isnan(grid)));
       hill = save_consts_struct.hill;
       
       if nan_sum ~= 0
           if hill == 5
               grid_nan(1) = grid_nan(1) + 1;
           elseif hill == 10
               grid_nan(2) = grid_nan(2) + 1;
           elseif hill == 50
               grid_nan(3) = grid_nan(3) + 1;
           elseif hill == 100
               grid_nan(4) = grid_nan(4) + 1;
           elseif hill == Inf
               grid_nan(5) = grid_nan(5) + 1;
           else
               error('Unexpected hill coef')
           end
       end
    end
end
% NaN only appears for hill-coef = 100!
%%
isequal(wave_vec_simple,wave_vec_fourier)

sum(wave_vec_simple ~= wave_vec_fourier)

sum(wave_vec_simple(wave_vec_simple ~= wave_vec_fourier))
%%
diag_count = 0;
for i = find(wave_vec_simple ~= wave_vec_fourier)
    if strcmp(orientation_vec{i},'Diagonal')
        diag_count = diag_count + 1;
    end
    
end
%% Matrix of all simulations and matrix mean for every hill coef
M_all = [hill_vec',wave_vec_simple',t_out_vec'];
%M_wave = M_all(:,2)
[C,ia,idx] = unique(M_all(:,1),'rows');
M_wave_frac = [C,accumarray(idx,M_all(:,2),[],@mean)];
%%
M_wave_frac2 = M_wave_frac;
M_wave_frac2(M_wave_frac2 == Inf) = 300;
figure
plot(M_wave_frac2(:,1),M_wave_frac2(:,2))
hold on
scatter(M_wave_frac2(:,1),M_wave_frac2(:,2), 'b')
ylim([0 1])
xlabel('Hill Coefficient')
ylabel('Fraction simulations wave')
xticks([10 50 100 300])
xticklabels({'10','50','100','Inf'})
%% plot including upper bound
M_wave_frac2 = M_wave_frac;
M_wave_frac2(M_wave_frac2 == Inf) = 300;
figure
plot(M_wave_frac2(:,1),M_wave_frac2(:,2))
hold on
plot(M_wave_frac2_rounded(:,1),M_wave_frac2_rounded(:,2))
scatter(M_wave_frac2(:,1),M_wave_frac2(:,2),[],[0    0.4470    0.7410])
scatter(M_wave_frac2_rounded(:,1),M_wave_frac2_rounded(:,2),[],[0.8500    0.3250    0.0980])
ylim([0 1])
xlabel('Hill Coefficient')
ylabel('Fraction simulations wave')
legend('Strict','Upper bound')
xticks([10 50 100 300])
xticklabels({'10','50','100','Inf'})


%% Plots about duration simulatoins
M_5 = M_all(M_all(:,1) == 5,:);
M_10 = M_all(M_all(:,1) == 10,:);
M_50 = M_all(M_all(:,1) == 50,:);
M_100 = M_all(M_all(:,1) == 100,:);
M_inf = M_all(M_all(:,1) == Inf,:);

M_vec = {M_5,M_10,M_50,M_100,M_inf};
M_vec_no_inf = {};

counter = 1;
for i = M_vec
    select_no_Inf = ~isinf(M_vec{counter}(:,3));
    M_vec_no_inf = [M_vec_no_inf, M_vec{counter}(select_no_Inf,:)];
    counter = counter + 1;
end

plot_duration = [ones(size(M_vec_no_inf{1}(:,3)));2* ones(size(M_vec_no_inf{2}(:,3))); 3* ones(size(M_vec_no_inf{3}(:,3))); 4* ones(size(M_vec_no_inf{4}(:,3))) ; 5 * ones(size(M_vec_no_inf{5}(:,3))) ];
%% boxplot
figure
boxplot([M_vec_no_inf{1}(:,3); M_vec_no_inf{2}(:,3); M_vec_no_inf{3}(:,3);M_vec_no_inf{4}(:,3);M_vec_no_inf{5}(:,3)],plot_duration)
set(gca,'XTickLabel',{'5','10','50','100','Inf'})
%% lineplot with std
mean_duration_vec = [mean(M_vec_no_inf{1}(:,3)) mean(M_vec_no_inf{2}(:,3)) mean(M_vec_no_inf{3}(:,3)) mean(M_vec_no_inf{4}(:,3)) mean(M_vec_no_inf{5}(:,3))];
std_duration_vec = [std(M_vec_no_inf{1}(:,3)) std(M_vec_no_inf{2}(:,3)) std(M_vec_no_inf{3}(:,3)) std(M_vec_no_inf{4}(:,3)) std(M_vec_no_inf{5}(:,3))];
%%
figure
l1 = plot([5 10 50 100 300],mean_duration_vec,'color', 'blue');
hold on
l2 = errorbar([5 10 50 100 300],mean_duration_vec,std_duration_vec,'color','blue');
xlabel('Hill coefficient')
ylabel('Duration simulation')
xticks([5 10 50 100 300])
xticklabels({'5','10','50','100','Inf'})
%% Plots t_onset (only waves)
M_wave = [hill_vec(wave_vec_simple == 1)',wave_vec_simple(wave_vec_simple == 1)',t_onset_vec(wave_vec_simple == 1)'];

M_5_wave =   M_wave(M_wave(:,1) == 5,:);
M_10_wave =  M_wave(M_wave(:,1) == 10,:);
M_50_wave =  M_wave(M_wave(:,1) == 50,:);
M_100_wave = M_wave(M_wave(:,1) == 100,:);
M_inf_wave = M_wave(M_wave(:,1) == Inf,:);

M_vec = {M_10_wave,M_50_wave,M_100_wave, M_inf_wave};

plot_duration = [ones(size(M_vec{1}(:,3)));2* ones(size(M_vec{2}(:,3))); 3* ones(size(M_vec{3}(:,3))); 4* ones(size(M_vec{4}(:,3)))];

%% boxplot
figure
boxplot([M_vec{1}(:,3); M_vec{2}(:,3); M_vec{3}(:,3); M_vec{4}(:,3)],plot_duration)
set(gca,'XTickLabel',{'10','50','100','Inf'})
xlabel('Hill coefficient')
ylabel('Time until wave')
%% I plots
indx = find(hill_vec == 10);

figure
for i = indx
    plot(1:length(I1_vec{i}),I1_vec{i})
    hold on
end
% --------------------------------------------------------------
%% Analyse random or wave initial situations
% note: cell states are called homogeneous when all cell states rounded to
% 6 decimals are identical. Less strict rounding sometimes says that
% simulations that stopped witout period before t_max are not homogeneous.
for hill_val = 1:12
    par_vec = []; % unique identifier for each parameter set
    C_on_vec = {};
    K_vec = {};
    period_vec = [];
    wave_strict_vec = []; % 1 if wave using strict requirements
    wave_rounded_vec = []; % 1 if wave using loose requirements (rounded p)
    t_out_vec = [];        
    t_onset_vec = [];
    homogeneous_oscillation_vec = []; % 1 if the pattern is a homogeneously oscillating grid
    homogeneous_end_vec = []; % 1 if final grid in simulation is homogeneous
    static_end_vec = []; % 1 if simulation ends in a static state (homogeneous or not)
    max_dif_hom_osc_vec = []; % max difference in cell state between time steps of period
    %% Obtain a vector containing paths to all parameter directories
    
    main_folder = sprintf('%s_%g','M:\tnw\bn\hy\Shared\Douwe\Recent Results\Hill_coef\Parameter_screening\Network_15\Hill',hill_val);
    %%
    dir_list = dir(main_folder);
    dir_list = dir_list([dir_list.isdir]);
    dir_list_all = dir_list(3:length(dir_list));
    dir_list_folder = {dir_list_all.folder};
    dir_list_name = {dir_list_all.name};
    dir_paths = {};
    
    for i = 1:length(dir_list_all)
        dir_paths = [dir_paths,sprintf('%s\\%s',dir_list_folder{i},dir_list_name{i})];
    end
    %% loop over all parameter folders and the simulations within each parameter
    % folder.
    for par_folder = 1:length(dir_paths)
        cd(dir_paths{par_folder})
        file_names = dir('*.mat');
        file_names = {file_names.name};
        
        rem_idx = strcmp(file_names,'parameters.mat');
        file_names(rem_idx) = [];
        %%
        for sim = 1:length(file_names)
            load(file_names{sim})
            
            % initialize values concerning homogeneous oscillations
            strict_static_end = 0;
            homogeneous_osc = 0; 
            max_homogeneous_osc = -1; %indicating no homogeneous oscillation.
            
            %% Check if end state is homogeneous
            std_end_g1 = std(cells_hist{length(cells_hist)}(:,1));
            std_end_g2 = std(cells_hist{length(cells_hist)}(:,2));
            
            std_hom = 0.02;
            
            if std_end_g1 < std_hom && std_end_g2 < std_hom
                homogeneous_end = 1;
            else
                homogeneous_end = 0;
            end
            
            %%
            % Old way:
            %{ 
             hom_end_g1 = all(round(cells_hist{length(cells_hist)}(:,1),6) == round(cells_hist{length(cells_hist)}(1,1),6));
             hom_end_g2 = all(round(cells_hist{length(cells_hist)}(:,2),6) == round(cells_hist{length(cells_hist)}(1,2),6));
                    
                    if hom_end_g1 == 1 && hom_end_g2 == 1
                        homogeneous_end = 1;
                    else
                        homogeneous_end = 0;
                    end
            %}
            %% Initialize wave values:
            wave_strict = 0;
            wave_rounded = 0;
            % If there is periodicity: check for wave & homogeneous
            % oscillations
            if t_onset ~= Inf
                orig_cells_hist = cells_hist;
                for t = 1:length(cells_hist)
                    cells_hist{t} = round(cells_hist{t},6);
                end
                
                wave_strict = travelling_wave_test(cells_hist, save_consts_struct.a0, period, t_out,distances);
                % Mostly using unrounded cells_states results in missing
                % waves, but sometimes you miss waves by rounding. So if
                % either of the two test below reveal wave, wave_rounded is
                % 1.
                wave_rounded1 = travelling_wave_test_rounded6(cells_hist, save_consts_struct.a0, period, t_out,distances);
                wave_rounded2 = travelling_wave_test_rounded6(orig_cells_hist, save_consts_struct.a0, period, t_out,distances);
                
                if wave_rounded1 ==1  || wave_rounded2 == 1
                    wave_rounded = 1;
                else
                    wave_rounded = 0;
                end
                
                % If in addition to a homogeneous end state, periodicity
                % is detected, there is a homogeneous oscillation. If the grid became
                % homogeneous because of a wave annhilation or a pattern
                % that died out, no periodicity would be detected. The end
                % state is homogeneous, so the other time steps of the
                % detected period must also be homogeneous.
                if period ~= Inf 
                    %%
                    state_g1_c1 = round(cells_hist{t_onset}(1,1),6);
                    state_g2_c1 = round(cells_hist{t_onset}(1,2),6);
                    
                    %% Check if homogeneous oscillations
                    std_onset_g1 = std(cells_hist{t_onset}(:,1));
                    std_onset_g2 = std(cells_hist{t_onset}(:,2));
                    
                    %Old:
                    %{
                    hom_g1 = sum(round(cells_hist{t_onset}(:,1),6) == state_g1_c1);
                    hom_g2 = sum(round(cells_hist{t_onset}(:,2),6) == state_g2_c1);
                    
                    if hom_g1 == size(cells_hist{1},1) && hom_g2 == size(cells_hist{1},1)
                    %}
                    if std_onset_g1 < std_hom && std_onset_g2 < std_hom && homogeneous_end == 1
                        %% Determine max difference between cell states
                        period_states = {};
                        for hom_state = 0:period - 1
                            period_states = [period_states, [round(cells_hist{t_onset + hom_state}(1,1),6) round(cells_hist{t_onset}(1,2),6)]];
                        end
                        
                        % pairwise differences
                        p2_val = 1;
                        max_dif = 0;
                        for p1 = 1:period - 1
                            p2_val = p2_val + 1;
                            for p2 = p2_val:period
                                dif = round(sum(abs(period_states{p1} - period_states{p2})),6);
                                if dif > max_dif
                                    max_dif = dif;
                                end
                            end
                        end
                        % If there is no difference in the four states of the
                        % homogeneous oscillations, there is no
                        % periodicity. then there is a static homogeneous
                        % grid, therefore, period = Inf. Rember, above
                        % selected for periodicity, if not found here, it
                        % means that only periodicity if not rounded to 6
                        % decimals.
                        if max_dif == 0
                            period = Inf;
                        else
                            homogeneous_osc = 1;
                            max_homogeneous_osc = max_dif;
                        end    
                    end
                end
            else
                if t_out < 10000
                        strict_static_end = 1;
                end
            end
            
            if strict_static_end == 1 || (homogeneous_osc == 0 && homogeneous_end == 1)
                static_end = 1;
            else
                static_end = 0;
            end
            
            % Store values
            period_vec = [period_vec,period];
            wave_strict_vec = [wave_strict_vec,wave_strict];
            wave_rounded_vec = [wave_rounded_vec,wave_rounded];
            t_out_vec = [t_out_vec,t_out];
            t_onset_vec = [t_onset_vec,t_onset];
            K_vec = [K_vec, save_consts_struct.K];
            C_on_vec = [C_on_vec, save_consts_struct.Con];
            par_vec = [par_vec,par_folder];
            max_dif_hom_osc_vec = [max_dif_hom_osc_vec,max_homogeneous_osc];
            homogeneous_oscillation_vec = [homogeneous_oscillation_vec,homogeneous_osc];
            homogeneous_end_vec = [homogeneous_end_vec,homogeneous_end];
            static_end_vec = [static_end_vec,static_end];
        end
        fprintf(sprintf('%g percent\n',par_folder/length(dir_paths)*100));
    end
    %% Save overview results above
    cd(main_folder)
    save(sprintf('new_Overview_Hill_%g_init_wave.mat',hill_val),'period_vec','wave_strict_vec','wave_rounded_vec','t_out_vec','t_onset_vec','K_vec','C_on_vec','par_vec','max_dif_hom_osc_vec','homogeneous_oscillation_vec','homogeneous_end_vec','dir_paths','static_end_vec')
    
    %% Prepare data for boxplot
    % Check if no other homogeneous oscillation periods etc
    true_wave_rounded = sum(wave_rounded_vec)/length(wave_rounded_vec);
    true_wave = sum(wave_strict_vec) / length(wave_strict_vec);
    dynamic_homogeneous = sum(homogeneous_oscillation_vec)/length(homogeneous_oscillation_vec);
    static_homogeneous = (sum(homogeneous_end_vec) - sum(homogeneous_oscillation_vec))/length(period_vec);
    static_end_state = sum(static_end_vec)/length(static_end_vec);
    infinite_dynamics = 1 - true_wave_rounded - dynamic_homogeneous - static_end_state;
    %% save
    cd ..
    save(sprintf('new_Boxplot_data_hill_%g',hill_val),'true_wave','dynamic_homogeneous','static_homogeneous','infinite_dynamics','true_wave_rounded', 'static_end_state');
end
%% ----------- Plotting boxplot
hill_vec =          [6 10];%12;
true_wave_vec =     zeros(1,2);%12);
dynamic_homogeneous_vec =  zeros(1,2);%12);
static_homogeneous_vec =  zeros(1,2);%12);
static_end_state_vec = zeros(1,2);
infinite_dynamics_vec =   zeros(1,2);%12);
%%
load_folder = 'M:\tnw\bn\hy\Shared\Douwe\Recent Results\Hill_coef\Network_15\Single_par_random_init\sims';%'W:\staff-groups\tnw\bn\hy\Shared\Douwe\Recent Results\Hill_coef\Parameter_screening\Network_19';
for i = [6 10]%1:12
    file_name = sprintf('new_Boxplot_data_hill_%g.mat',i);
    load(fullfile(load_folder, file_name))
    
    true_wave_vec(i) =              true_wave_rounded;
    dynamic_homogeneous_vec(i) =    dynamic_homogeneous;
    static_homogeneous_vec(i) =     static_homogeneous;
    infinite_dynamics_vec(i) =      infinite_dynamics;
    static_end_state_vec(i) =       static_end_state;
end
%%
M_all = [true_wave_vec' infinite_dynamics_vec' dynamic_homogeneous_vec' static_end_state_vec'];

figure
bar(M_all,'stacked')
legend('Pure wave','Infinite Dynamics', 'Homogeneous oscillations','Static','location','northeastoutside')
xlabel('Hill coefficient')
ylabel('Fraction of tried parameters')
%% Plot maximal difference between homogeneous oscillation grids
max_dif_vec = {};
mean_max_dif_vec = zeros(1,12);
std_max_dif_vec = zeros(1,12);

mean_t_out_vec = zeros(1,12);
std_t_out_vec = zeros(1,12);

for i = 1:12
    load_path = sprintf('W:\\staff-groups\\tnw\\bn\\hy\\Shared\\Douwe\\Recent Results\\Hill_coef\\Parameter_screening\\Network_19\\Hill_%g\\Overview_Hill_%g_init_wave',i,i);
    load(load_path)
    
    max_dif_hom_osc_vec(max_dif_hom_osc_vec == -1) = [];
   
    max_dif_vec = [max_dif_vec, max_dif_hom_osc_vec];
    mean_max_dif_vec(i) = mean(max_dif_hom_osc_vec);
    std_max_dif_vec(i) = std(max_dif_hom_osc_vec);
    
    % calculate how long it takes before patterns form, so not taking into 
    % account simulations that did not form a pattern
    t_out_vec( find(t_out_vec == Inf)) = NaN;
    mean_t_out_vec(i) = nanmean(t_out_vec);
    std_t_out_vec(i) = nanstd(t_out_vec);
end


%% Max dif homogeneous oscillations
figure
plot(1:12,mean_max_dif_vec)
hold on
errorbar(1:12,mean_max_dif_vec,std_max_dif_vec,'color','blue')
xlim([0 12.5])
ylim([0 2])
xlabel('Hill coefficient')
ylabel('Max distance homogeneously oscillating grids')
%% mean duration
figure
plot(1:12,mean_t_out_vec)
hold on
errorbar(1:12,mean_t_out_vec,std_t_out_vec,'color','blue')
xlim([0 12.5])
%ylim([0 2])
xlabel('Hill coefficient')
ylabel('Mean duration simulations')
%% Checking examples
inf_dym = homogeneous_end_vec == 0 & wave_rounded_vec == 0;
indx_inf_dym = find(inf_dym == 1);
inf_dym_sims = dir_paths(indx_inf_dym);

stat_homs = homogeneous_end_vec == 1 & homogeneous_oscillation_vec == 0;
indx_stat_homs = find(stat_homs == 1);
stat_hom_sims = dir_paths(indx_stat_homs);




