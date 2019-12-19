function update_cell_figure_external_temp(h_cells, h_borders, cells, t, disp_mol, pos)  
    % plot title
    title(sprintf('t=%d', t), 'FontSize', 20);
    
    % Update cell colors
    %cells = ones(N,1);
    N = size(cells, 1);
    c_all = zeros(N, 3);
    c_all(:, 2) = cells; % green channel
    set(h_cells, 'cdata', c_all);

    % Update cell positions
    set(h_cells, 'xdata', pos(:, 1), 'ydata', pos(:, 2) );
    %set(h_borders, 'xdata', pos(:, 1), 'ydata', pos(:, 2) );
end