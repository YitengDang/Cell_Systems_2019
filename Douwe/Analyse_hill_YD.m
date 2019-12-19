clear all 
close all
clc
set(0, 'defaulttextinterpreter', 'tex');
%% Load analyzed data

% TW formation (from random initial conditions)
%{
load_folder = 'M:\tnw\bn\hy\Shared\Douwe\Recent Results\Hill_coef\Network_15\Single_par_random_init\sims';
hill_all = [1:12 Inf]; 
label = 'TW_formation';
%}

% TW propagation
%
load_folder = 'M:\tnw\bn\hy\Shared\Douwe\Recent Results\Hill_coef\Parameter_screening\Network_15';%'W:\staff-groups\tnw\bn\hy\Shared\Douwe\Recent Results\Hill_coef\Parameter_screening\Network_19';
hill_all = [1:12]; 
label = 'TW_propagation';
%}
% If from remote location:
%load_folder = strrep(load_folder, 'M:\', 'W:\staff-groups\');
true_wave_vec = zeros(numel(hill_all), 1);
periodic_vec = zeros(numel(hill_all), 1);
static_vec = zeros(numel(hill_all), 1);
infinite_dynamics_vec = zeros(numel(hill_all), 1);
frac_dynamic = zeros(numel(hill_all), 1);
for i = 1:numel(hill_all)
    file_name = sprintf('Boxplot_data_hill_%g_YD.mat', hill_all(i) );
    load(fullfile(load_folder, file_name))
    
    % Yiteng's classification
    true_wave_vec(i) =              true_wave_rounded;
    periodic_vec(i) =               periodic;
    static_vec(i) =                 static;
    infinite_dynamics_vec(i) =      infinite_dynamics;
    
    % Binary classification
    frac_dynamic(i) = true_wave_rounded+periodic+infinite_dynamics;
    % Douwe's classification
    %dynamic_homogeneous_vec(i) =    dynamic_homogeneous;
    %static_homogeneous_vec(i) =     static_homogeneous;
    %infinite_dynamics_vec(i) =      infinite_dynamics;
end
%% Plot analyzed results: binary classification

% TW formation
%{
x_data = hill_all; x_data(x_data==Inf) = hill_all(end-1)+2;
y_data = frac_dynamic';
% reverse data
x_data = hill_all(hill_all~=Inf); x_data = [-1 x_data];
xticklabels = sprintfc('%d', flip(x_data(2:end))); xticklabels = ['\infty' xticklabels];
y_data = flip(frac_dynamic');
%}

% TW propagation
x_data = hill_all;
y_data = frac_dynamic';
xticklabels = sprintfc('%d', x_data); 
y_data = flip(frac_dynamic');
%}

h=figure;
plot(x_data, y_data, 'bo-', 'LineWidth', 1.5);
xlabel('Hill coefficient')
set(gca, 'XTick', x_data, 'XTickLabels',xticklabels);
ylabel('Fraction of simulations')
set(gca, 'FontSize', 32);
box on
xlim([x_data(1)-1/2 x_data(end)+1/2])
ylim([0 1]);

% Save figure
qsave = 0;
folder = 'H:\My Documents\Multicellular automaton\figures\two_signals\trav_wave_vs_Hill';
fname_str = sprintf('%s_fraction_dynamic_patterns_simulations_YD_flipped_size_10_8', label); %_fontsize32_fig_size10_8';
fname = fullfile(folder, fname_str);
save_figure(h, 10, 8, fname, '.pdf', qsave);

% save figure data
if qsave
    save(fname, 'x_data', 'y_data');
end
%% Plot analyzed results: Fraction TWs

% TW formation
%x_data = [hill_all(1:end-1) hill_all(end-1)+2];
%y_data = true_wave_vec'; 
% TW propagation
x_data = 1:numel(hill_all);
y_data = true_wave_vec';
%x_data = [-1 1:numel(hill_all)];
%y_data = [true_wave_vec' 1];
%xticklabels = sprintfc('%d', x_data); %{'1','2','3','4','5','6','7','8','9','10','11','12','\infty'}

% reverse data
%xticklabels = cell(numel(hill_all)+1, 1);
%xticklabels{1} = '\infty';
%xticklabels(2:end) = sprintfc('%d', flip(hill_all)); 
xticklabels = sprintfc('%d', flip(hill_all)); 

y_data = flip(y_data);

h=figure;
hold on
plot(x_data, y_data, 'b--', 'LineWidth', 1.5);
scatter(x_data, y_data, 100, 'b', 'filled');

xlabel('Hill coefficient')
set(gca, 'XTick', x_data, 'XTickLabels', xticklabels);
ylabel('Fraction of simulations')
set(gca, 'FontSize', 32);
box on
xlim([1 12]);
%ylim([0 0.3]);
ylim([0 1]);

% Save figure
qsave = 1;
folder = 'H:\My Documents\Multicellular automaton\figures\two_signals\trav_wave_vs_Hill';
fname_str = sprintf('%s_fraction_TWs_simulations_YD_reversed_size_10_8_v3', label); %_fontsize32_fig_size10_8';
fname = fullfile(folder, fname_str);
save_figure(h, 10, 8, fname, '.pdf', qsave);

% save figure data
save_data = 0;
if qsave && save_data
    save(fname, 'x_data', 'y_data');
end
%% Plot analyzed results: 4-way classification
%M_all = [true_wave_vec' infinite_dynamics_vec' dynamic_homogeneous_vec' static_homogeneous_vec'];
%M_all = [true_wave_vec' periodic_vec' static_vec' infinite_dynamics_vec'];
M_all = [true_wave_vec periodic_vec static_vec infinite_dynamics_vec];

x_data = hill_all;
bar_data = M_all;

h=figure;
bar(bar_data,'stacked')
%legend('Pure wave','Infinite Dynamics', 'Homogeneous oscillations','Static Homogeneous','location','northeastoutside')
legend('Travelling wave','Periodic', 'Static','Infinite Dynamics','location','nw'); %'northeastoutside')
xlabel('Hill coefficient')
xticklabels({'1','2','3','4','5','6','7','8','9','10','11','12','\infty'})
ylabel('Fraction of simulations')
set(gca, 'FontSize', 32);
box on

% save figure
qsave = 0;
%folder = 'W:\staff-groups\tnw\bn\hy\Shared\Douwe\Recent Results\Hill_coef\Parameter_screening\Network_15';
%folder = 'M:\tnw\bn\hy\Shared\Douwe\Recent Results\Hill_coef\Parameter_screening\Network_15';
%folder = 'M:\tnw\bn\hy\Shared\Douwe\Recent Results\Hill_coef\Network_15\Single_par_random_init\sims';
folder = 'H:\My Documents\Multicellular automaton\figures\two_signals\trav_wave_vs_Hill';

fname_str = sprintf('%s_barplot_endstate_simulations_YD', label); %_fontsize32_fig_size10_8';
fname = fullfile(folder, fname_str);
save_figure(h, 12, 8, fname, '.pdf', qsave);

% Save data 
if qsave
    save(fname, 'x_data', 'bar_data');
end
%% Analyse random or wave initial situations
% note: cell states are called homogeneous when all cell states rounded to
% 6 decimals are identical. Less strict rounding sometimes says that
% simulations that stopped witout period before t_max are not homogeneous.
hill_all = [1:12];
for hill_val = hill_all
    %%
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
    max_dif_hom_osc_vec = []; % max difference in cell state between time steps of period
    %% Obtain a vector containing paths to all parameter directories
    %main_folder = sprintf('%s_%g','M:\tnw\bn\hy\Shared\Douwe\Recent Results\Hill_coef\Parameter_screening\Network_15\Hill',hill_val);
    %main_folder = sprintf('%s_%g','W:\staff-groups\tnw\bn\hy\Shared\Douwe\Recent Results\Hill_coef\Parameter_screening\Network_15\Hill',hill_val);
    
    main_folder = sprintf('%s_%g', 'M:\tnw\bn\hy\Shared\Douwe\Recent Results\Hill_coef\Network_15\Single_par_random_init\sims\Hill', hill_val);
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
    % folder
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
        end
        fprintf(sprintf('%g percent\n',par_folder/length(dir_paths)*100));
    end
    %}
    %% Save overview results above
    %
    cd(main_folder)
    save(sprintf('Overview_Hill_%g_init_wave_YD.mat',hill_val),...
        'period_vec','wave_strict_vec','wave_rounded_vec','t_out_vec','t_onset_vec',...
        'K_vec','C_on_vec','par_vec','max_dif_hom_osc_vec',...
        'homogeneous_oscillation_vec','homogeneous_end_vec','dir_paths');
    %}
    %% Load overview results
    %folder = fullfile('W:\staff-groups\tnw\bn\hy\Shared\Douwe\Recent Results\Hill_coef\Parameter_screening\Network_15',...
    %    sprintf('Hill_%d', hill_val) );
    %folder = fullfile('W:\staff-groups\tnw\bn\hy\Shared\Douwe\Recent Results\Hill_coef\Network_15\Single_par_random_init\sims',...
    %    sprintf('Hill_%d', hill_val) );
    
    %folder = fullfile('M:\tnw\bn\hy\Shared\Douwe\Recent Results\Hill_coef\Parameter_screening\Network_15',...
    %    sprintf('Hill_%d', hill_val) );
    %folder = fullfile('M:\tnw\bn\hy\Shared\Douwe\Recent Results\Hill_coef\Network_15\Single_par_random_init\sims',...
    %    sprintf('Hill_%d', hill_val) );
    folder = main_folder;
    
    load(fullfile(folder, sprintf('Overview_Hill_%g_init_wave.mat',hill_val)),...
        'period_vec','wave_strict_vec','wave_rounded_vec','t_out_vec','t_onset_vec',...
        'K_vec','C_on_vec','par_vec','max_dif_hom_osc_vec',...
        'homogeneous_oscillation_vec','homogeneous_end_vec','dir_paths');
    
    %% Prepare data for boxplot
    % Classification YD
    t_max = 10000; % from simulation data
    true_wave_rounded = sum(wave_rounded_vec)/length(wave_rounded_vec);
    true_wave = sum(wave_strict_vec) / length(wave_strict_vec);
    static = sum(period_vec==Inf & t_out_vec<t_max)/length(period_vec);
    periodic = sum(period_vec<Inf & ~wave_rounded_vec)/length(period_vec);
    infinite_dynamics = 1 - true_wave_rounded - static - periodic;
    %%
    % Classification DG
    % Check if no other homogeneous oscillation periods etc
    %{
    true_wave_rounded = sum(wave_rounded_vec)/length(wave_rounded_vec);
    true_wave = sum(wave_strict_vec) / length(wave_strict_vec);
    dynamic_homogeneous = sum(homogeneous_oscillation_vec)/length(homogeneous_oscillation_vec);
    static_homogeneous = (sum(homogeneous_end_vec) - sum(homogeneous_oscillation_vec))/length(period_vec);
    infinite_dynamics = 1 - true_wave_rounded - dynamic_homogeneous - static_homogeneous;
    %}
    % save
    %cd ..
    %save(sprintf('Boxplot_data_hill_%g_YD',hill_val),...
    %    'true_wave','true_wave_rounded','dynamic_homogeneous','static_homogeneous','infinite_dynamics');
    %folder = fullfile('W:\staff-groups\tnw\bn\hy\Shared\Douwe\Recent Results\Hill_coef\Parameter_screening\Network_15');
    %folder = fullfile('M:\tnw\bn\hy\Shared\Douwe\Recent Results\Hill_coef\Parameter_screening\Network_15');
    %folder = 'W:\staff-groups\tnw\bn\hy\Shared\Douwe\Recent Results\Hill_coef\Network_15\Single_par_random_init\sims';
    %folder = 'M:\tnw\bn\hy\Shared\Douwe\Recent Results\Hill_coef\Network_15\Single_par_random_init\sims';
    cd(main_folder)
    cd ..
    save( sprintf('Boxplot_data_hill_%g_YD', hill_val),...
        'true_wave','true_wave_rounded','static','periodic', 'infinite_dynamics');
end
%% ----------- Plotting boxplot
end_val = numel(hill_all);
hill_vec = 1:end_val;
%{|
true_wave_vec =     zeros(1,end_val);
dynamic_homogeneous_vec =  zeros(1,end_val);
static_homogeneous_vec =  zeros(1,end_val);
infinite_dynamics_vec =   zeros(1,end_val);
%}
true_wave_vec = zeros(1,end_val);
periodic_vec = zeros(1,end_val);
static_vec = zeros(1,end_val);
infinite_dynamics_vec = zeros(1,end_val);

%{
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
%}