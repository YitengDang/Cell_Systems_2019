

% load original initial state and other required variables
folder = 'D:\Documenten\Studie\Master\Internship Delft\selected patterns\Snapshots';
fname_str = 'wave_horizontal_downward';
load(fullfile(folder, fname_str), 'cells_hist', 'distances', 'save_consts_struct');
a0 = save_consts_struct.a0;
Con = save_consts_struct.Con;
Coff = save_consts_struct.Coff;
K = save_consts_struct.K;
gridsize = sqrt(save_consts_struct.N);
M_int = save_consts_struct.M_int;
lambda = save_consts_struct.lambda;

% set random seed
seed = rng('shuffle'); 
%% Alter original initial state and run simulation
%N = size(distances,1);
N = 1;
orig_init = cells_hist{1};
init_state = orig_init;

% indicate whether simulations need to be saved in 'app format'
save_app = 0;

%% Change each cell to every possible state, plane waves
changed_cell_vec = [];
state_to_vec = [];
state_from_vec = [];
t_out_vec = [];
t_wave_vec = [];
orientation_vec = {};
direction_vec = {};
nr_waves_vec = [];
nr_bands_vec = [];
bended_vec = [];
diag_wrap_vec = {};

for cell_nr = 1:N
    for state = 1:4 % state 1 = [0 0], 2 = [1 0], 3 = [0 1], 4 = [1 1]
        init_state = orig_init;
        if ~isequal(orig_init(cell_nr,:),translate_states(state))
            init_state(cell_nr,:) = translate_states(state);
            [cells_hist,synchrony_hist1,synchrony_hist2,dist,t_wave,t_out]...
                = time_evolution_2_signals_func(a0,Con,Coff,K,gridsize,init_state,M_int,lambda);
            if t_wave ~= -1
                [orientation,bands_in_wave,number_of_waves,diag_wrap,...
                    band_vec,wave_state,bended] = determine_wave_properties(cells_hist,t_wave);
                [~,~,~,~,band_vec2, ~,~] = determine_wave_properties(...
                    cells_hist,t_wave + 1,wave_state);
                if strcmp(orientation,'Horizontal') || strcmp(orientation,'Vertical')
                    direction = Determine_travel_direction(band_vec,...
                        band_vec2,orientation,bands_in_wave);
                elseif strcmp(orientation,'Diagonal')
                    direction_hor = Determine_travel_direction(...
                        band_vec{1},band_vec2{1},'Horizontal',sum(band_vec{1}));
                    direction_ver = Determine_travel_direction(...
                        band_vec{2},band_vec2{2},'Vertical',1);
                    direction = sprintf('%s & %s',direction_hor,direction_ver);
                else
                    direction = NaN;
                end
            else % in case no wave is detected
                orientation = -1;
                bands_in_wave = -1;
                number_of_waves = -1;
                diag_wrap = -1;
                band_vec = -1;
                wave_state = -1;
                bended = -1;
                band_vec2 = -1;
                direction = -1;
            end
            % Store values
            changed_cell_vec = [changed_cell_vec,cell_nr];
            state_to_vec = [state_to_vec, state];
            state_from_vec = [state_from_vec,translate_states(orig_init(cell_nr,:))];
            t_out_vec = [t_out_vec,t_out];
            t_wave_vec = [t_wave_vec,t_wave];
            orientation_vec = [orientation_vec,orientation];
            direction_vec = [direction_vec, direction];
            nr_waves_vec = [nr_waves_vec,number_of_waves];
            nr_bands_vec = [nr_bands_vec,bands_in_wave];
            bended_vec = [bended_vec,bended];
            diag_wrap_vec = [diag_wrap_vec,diag_wrap];
            
            if save_app == 1
                save_app_format(folder,fname_str,cells_hist,cell_nr,state,...
                    synchrony_hist1,synchrony_hist2)
            end
        end
    end
end

% save results VOEG ALLE VECS NOG TOE
%save(sprintf('%s_overview_robustness_single_change.mat',fname_str),'changed_cell_vec','state_to_vec','state_from_vec','t_out_vec','t_wave_vec');
%% plot time state heatmap
plot_time_state_heatmap(cells_hist,12)