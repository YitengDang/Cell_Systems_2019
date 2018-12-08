function [status, cells, fname_ini_state] = manual_input_state(...
    signal_count, folder, N)
    % used in simulate_cells_... functions when manual input state is called
    disp('Manually input initial state...');

    %folder = app.InputStateSaveFolder; %fullfile('..', '..', 'data', 'input_states');
    [file, path] = uigetfile({'*.xls'; '*.xlsx'; '*.*'},...
        'Load initial lattice configuration', folder );
    [~, fname_ini_state, ~] = fileparts(file);
    
    try
        cells_in = cell(2, 1);
        cells_in{1} = xlsread(fullfile(path, file), 1);
        cells_in{2} = xlsread(fullfile(path, file), 2);
    catch
        try 
            cells_in = cell(2, 1);
            cells_in{1} = xlsread(fullfile(path, file), 1);
        catch 
            error('Wrong input initial state, operation aborted');
            status = 0;
            cells = [];
            return
        end
    end

    [status, cells] = InputCellStateParser(cells_in, signal_count, N);
    if ~status
        error('Wrong input data format');
        return
    end
end
        
function [status, cells_parsed] = InputCellStateParser(cells, signal_count, N)
    % Checks whether input cells state (loaded from xlsx file) is in the right format 
    % If so, parse cell states
    %[gz, ~] = get_input_vars(app);
    %N = gz^2;
    gz = sqrt(N);
    status = 1; % Default
    if signal_count==1
        % one signal
        cells = cells{1};
        format1 = (size(cells, 1)==N && size(cells,2)==1);
        format2 = (size(cells, 1)==gz && size(cells,2)==gz);
        if format1 % N x 1 arrangement
            cells_parsed = cells;
        elseif format2 % gz x gz arrangement
            cells_parsed = cells(:);
        else
            cells_parsed = [];
            status = 0; % cells not good
        end
    elseif signal_count==2
        % two signals
        format1 = (size(cells{1}, 1)==N && size(cells{1},2)==2);
        format2 = ( all(size(cells{1})==[gz gz]) && all(size(cells{2})==[gz gz]) );
        if format1 % N x 2 arrangement
            cells_parsed = cells{1};
        elseif format2 % two sheets, gz x gz arrangement
            cells_parsed = zeros(N, 2);
            cells_parsed(:, 1) = reshape(cells{1}, N, 1);
            cells_parsed(:, 2) = reshape(cells{2}, N, 1);
        else
            cells_parsed = [];
            status = 0; % cells not good
        end
    end
end