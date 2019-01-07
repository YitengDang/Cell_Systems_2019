function[orientation,bands_in_wave,number_of_waves,diag_vec,band_vec,...
    wave_state,bended,wave_states] = determine_wave_properties(...
    cells_hist,t_wave,wave_state)
% Determine whether the wave is horizontal, vertical , or diaganol oriented
% on the grid. For horizontal and vertical orientation the number of waves
% and the number of bands within a wave is also determined. For diagonal
% oriented waves this is not determined. Instead Diagonal waves are 
% described in terms of how many times the wave wraps around the x and y 
% axis.
% This function should only be used when there is actually a wave present, 
% otherwise the output makes no sense.  
% This function also determines whether waves are bended ('met een knik')
% or not. If waves are bended, some other outputs may be unrealisable in
% some cases. For vertical and horizontal orientation, the number of waves, 
% and the number of bands may be incorect, however in the examples so far
% everything is fine. For diagonal orientation the travel direction may not
% be determined correctly if there are many bends.

% The wave state can be given as an input to this function to select a
% specific waveband colour. Usefull when calculating band_vec at t + 1.

% clarification inputs:
% cells_hist -> cell states of simulation
% t_wave -> time when wave was formed (determine using travelling_wave_test)
% wave_state -> see explanatino above

% clarification outputs:
% orientation -> how wave are oriented on grid: horizontal, vertical,diagonal
% bands_in_wave -> number of bands of specific colour in a single wave 
% diag_vec -> diagonal wave wraps around axis in vector [x_wraps y_wraps]
% band_vec -> vector with 1 at position band
% wave_state -> a cell state of the wave
% bended -> if wave is bended. 1 if yes, 0 if not

if nargin < 3
    wave_state = -1; % Indicate that wave_state can be freely chosen
end
%%
cells_hist_translated = zeros(size(cells_hist{1},1),1);
gridsize = sqrt(size(cells_hist_translated,1));

for i = 1:size(cells_hist_translated,1)
    cells_hist_translated(i) = translate_states(cells_hist{t_wave + 1}(i,:)); % +1 because at index 1 is t = 0
end

%% Identify background state
background_state = mode(cells_hist_translated);
cells_hist_translated(cells_hist_translated == background_state) = 0;
%% Pick one wave colour
wave_states = unique(cells_hist_translated(cells_hist_translated >0));
if wave_state == -1
    wave_state = datasample(wave_states,1);
end
wave_state_cell_idx = find(cells_hist_translated == wave_state);

%% Each row in the two matrices below holds the cell indices of either a 
% row or a column of the original grid. A row in M_rows corresponds to a row
% in the original grid and a row in M_columns corresponds to a column in the
% original grid. There is a second M_column, because waves can be oriented
% with 'left zigzag' instead of 'right zigzag'.
M_rows = zeros(gridsize,gridsize);
M_columns1 = zeros(gridsize,gridsize); % right zigzag
M_columns2 = zeros(gridsize,gridsize); % left zigzag

% Fill in cell indices for left zigzag
M_columns2(:,1) = 1:gridsize:gridsize*gridsize-1;
last_col_right_zig = 1+gridsize*(gridsize-1):1:1+gridsize*(gridsize-1)+(gridsize-1);
first_col_right_zig = 1:1:1+(gridsize-1);
for i = 2:gridsize
    if mod(i,2) ==0
        M_columns2(1,i) = last_col_right_zig(i);
    else
        M_columns2(1,i) = first_col_right_zig(i);
    end
end

for i = 2:gridsize
    zigzag = -1;
    for j = 2:gridsize
        if zigzag == 1
            cell_idx = M_columns2(i,j-1) + (gridsize + 1);
        elseif zigzag == -1
            cell_idx = M_columns2(i,j-1) - (gridsize - 1);
        end
        
        if cell_idx < 0
            cell_idx = last_col_right_zig(j);
        end
        M_columns2(i,j) = cell_idx;
        zigzag = zigzag * -1;
    end
end

% Entries in row_req or col_req are 1 if that row/column of the original 
% grid have a wave_state cell. Row_band and col_band only have a 1 if at
% the number of cells in that row/columnis in wave_state have crossed a
% treshold. Hereby bends can be identified.
row_req = zeros(1,gridsize);
col_req1 = zeros(1,gridsize);
col_req2 = zeros(1,gridsize);
row_band = zeros(1,gridsize);
col_band1 = zeros(1,gridsize);
col_band2 = zeros(1,gridsize);

for i = 1:gridsize
    % Placing the correct grid cell indices in each matrix row
    M_rows(i,:) = i:gridsize:i+(gridsize - 1)*gridsize;
    M_columns1(i,:) = 1+gridsize*(i-1):1:1+gridsize*(i-1)+(gridsize-1);
    
    % Checking if row/column is required to cover the wave state cells
    row_intersect = intersect(M_rows(i,:),wave_state_cell_idx);
    col_intersect1 = intersect(M_columns1(i,:),wave_state_cell_idx);
    col_intersect2 = intersect(M_columns2(i,:),wave_state_cell_idx);
    
    % Rows:
    if ~isempty(row_intersect)
        row_req(i) = 1;
        if length(row_intersect) >= 0.5 * gridsize
            row_band(i) = 1;
        end
    end
    
    % Columns:
    if ~isempty(col_intersect1)
        col_req1(i) = 1;
        if length(col_intersect1) > round(0.5 * (gridsize -1) * 0.6)
            col_band1(i) = 1;
        end
    end
    
    if ~isempty(col_intersect2)
        col_req2(i) = 1;
        if length(col_intersect2) > round(0.5 * (gridsize -1) * 0.6)
            col_band2(i) = 1;
        end
    end
end

% Decide if right or left zigzag columns best capture the wave
if sum(col_req1) < sum(col_req2)
    col_req = col_req1;
    col_band = col_band1;
    M_columns = M_columns1;
else
    col_req = col_req2;
    col_band = col_band2;
    M_columns = M_columns2;
end

%% Determine orientation 
sum_row_req = sum(row_req);
sum_col_req = sum(col_req);

if sum_row_req == sum_col_req
    orientation = 'Diagonal';
elseif sum_row_req < sum_col_req
    orientation = 'Horizontal';
else
    orientation = 'Vertical';
end
%% Determine the number of bands for each wave colour (assumed to be the
% same for each colour). Multiple bands can be caused by multiple waves, a
% single wave with multiple bands, or multiple waves with multiple bands.
% Also, determine whether the wave is bended or not.

diag_vec = NaN; 

if strcmp(orientation,'Horizontal')
    col1 = M_columns(1,:);
    col2 = M_columns(2,:);
    bands1 = intersect(col1,wave_state_cell_idx);
    bands2 = intersect(col2,wave_state_cell_idx);
    
    band1_vec = zeros(length(col1),1);
    band2_vec = zeros(length(col2),1);
    for i = 1:length(col1)
        band1_vec(i) = ismember(col1(i),bands1);
        band2_vec(i) = ismember(col2(i),bands2);
    end
   
    if sum(band1_vec) > sum(band2_vec)
        band_vec = band1_vec;
    else
        band_vec = band2_vec;
    end
    
    [bands_in_wave,number_of_waves] = determine_band_width(row_band);

    if row_req == row_band
        bended = 0;
    else
        bended = 1;
    end
    
elseif strcmp(orientation,'Vertical')
    row1 = M_rows(1,:);
    row2 = M_rows(2,:);
    bands1 = intersect(row1,wave_state_cell_idx);
    bands2 = intersect(row2,wave_state_cell_idx);
    
    band1_vec = zeros(length(row1),1);
    band2_vec = zeros(length(row2),1);
    for i = 1:length(row1)
        band1_vec(i) = ismember(row1(i),bands1);
        band2_vec(i) = ismember(row2(i),bands2);
    end
   
    if sum(band1_vec) > sum(band2_vec)
        band_vec = band1_vec;
    else
        band_vec = band2_vec;
    end
    
    [bands_in_wave,number_of_waves] = determine_band_width(band_vec);
   
    if col_req == col_band
        bended = 0;
    else
        bended = 1;
    end
    
elseif strcmp(orientation,'Diagonal')
    bands_in_wave = NaN;
    number_of_waves = NaN;
    y_wrap = length(intersect(M_rows(1,:),wave_state_cell_idx));
    x_wrap = length(intersect(M_columns(1,:),wave_state_cell_idx));
    diag_vec = [x_wrap y_wrap];
    bended = NaN;
    
    % horizontal band vec
    col1 = M_columns(1,:);
    bands_hor = intersect(col1,wave_state_cell_idx);
    band_vec_hor = zeros(length(col1),1);
    for i = 1:length(col1)
        band_vec_hor(i) = ismember(col1(i),bands_hor);
    end
    % vertical band vec
    row1 = M_rows(1,:);
    bands_ver = intersect(row1,wave_state_cell_idx);
    band_vec_ver = zeros(length(row1),1);
    for i = 1:length(row1)
        band_vec_ver(i) = ismember(row1(i),bands_ver);
    end
    
    band_vec = {band_vec_hor,band_vec_ver};
else
    error('Undetermined wave orientation.')
end



end