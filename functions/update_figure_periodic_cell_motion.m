function update_figure_periodic_cell_motion(h_cells, h_borders, cells, t, disp_mol, showI, a0, dist, pos)
% For moving cells, replot the lattice at each time step
    if nargin < 4
        disp_mol = 1;
    elseif nargin<5
        showI=0;
    end
    %% 
    if disp_mol==12
        % no signal -> white
        % signal 1 -> yellow
        % signal 2 -> blue 
        % signals 1&2 -> black

        % plot title
        if showI %whether or not to display I (slow down simulation)
            this_p = mean(cells, 1);
            this_I = [moranI(cells(:, 1), a0*dist) moranI(cells(:, 2), a0*dist)];
            title(sprintf('t=%d, p1 = %.2f, p2 = %.2f, I1 = %.2f, I2 = %.2f', t,...
                this_p(1), this_p(2), this_I(1), this_I(2)), ...
                'FontSize', 20);
        else
            this_p = mean(cells, 1);
            title(sprintf('t=%d, p1 = %.2f, p2 = %.2f', t, this_p(1), this_p(2)), ...
                'FontSize', 20);
        end

        % --update cells--
        
        % Update cell states
        clrs = 1-cells;
        c_all = zeros(size(cells, 1), 3); 
        c_all(:, 3) = clrs(:, 1); % signal 1 present -> Turn on blue channel
        c_all(:, 2) = clrs(:, 2); % signal 2 present -> Turn on green channel
        c_all(:, 1) = clrs(:, 2); % signal 2 present -> Turn on red channel
        set(h_cells, 'cdata', c_all);
        
        % Update cell positions
        set(h_cells, 'xdata', pos(:, 1), 'ydata', pos(:, 2) );
        set(h_borders, 'xdata', pos(:, 1), 'ydata', pos(:, 2) );
    else
        cells = cells(:, disp_mol);
        
        % Update title
        i_str = '';
        if showI 
            this_p = mean(cells);
            this_I = moranI(cells, a0*dist);
            title(sprintf('t=%d, p%s = %.2f, I%s = %.2f', t, i_str, this_p, i_str, this_I), ...
                'FontSize', 20);
        else
            this_p = mean(cells);
            title(sprintf('t=%d, p%s = %.2f', t, i_str, this_p), ...
                'FontSize', 20);
        end
        
        % Update cell colors
        %cells = ones(N,1);
        c_all = repmat(1-cells, 1, 3);
        set(h_cells, 'cdata', c_all);
        
        % Update cell positions
        set(h_cells, 'xdata', pos(:, 1), 'ydata', pos(:, 2) );
        set(h_borders, 'xdata', pos(:, 1), 'ydata', pos(:, 2) );
    end
end