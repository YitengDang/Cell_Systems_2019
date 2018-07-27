function update_cell_figure_continuum(plothandle, dist, a0, cells, t, disp_mol, showI)

    ax = gca;
    set(ax, 'FontSize', 20);
    if disp_mol==12 % plot both colors
        % no signal -> white
        % signal 1 -> yellow
        % signal 2 -> blue 
        % signals 1&2 -> black

        % plot title
        
        title(ax, sprintf('$$t=%d$$', t));
        %
        if showI %whether or not to display I (slow down simulation)
            %title(ax, sprintf('t=%d, p = %.2f', t, mean(cells)), ...
            %    'FontSize', 20);
            this_p = mean(cells, 1);
            this_I = [moranI(cells(:, 1), a0*dist) moranI(cells(:, 2), a0*dist)];
            title(ax, sprintf('$$t=%d, p_1 = %.2f, p_2 = %.2f, I_1 = %.2f, I_2 = %.2f$$', t,...
                this_p(1), this_p(2), this_I(1), this_I(2)), ...
                'FontSize', 20);
        else
            this_p = mean(cells, 1);
            title(ax, sprintf('$$t=%d, p_1 = %.2f, p_2 = %.2f$$', t, this_p(1), this_p(2)), ...
                'FontSize', 20);
        end
        %}
        % --update cells--
        clrs = 1-cells;
        c_all = zeros(size(cells, 1), 3); 
        c_all(:, 3) = clrs(:, 1); % signal 1 present -> Turn on blue channel
        c_all(:, 2) = clrs(:, 2); % signal 2 present -> Turn on green channel
        c_all(:, 1) = clrs(:, 2); % signal 2 present -> Turn on red channel
    else
        % plot title
        title(ax, sprintf('$$t=%d$$', t));
        %
        % show molecule number (not optimized)
        i_str = ''; 
        if ~isempty(disp_mol)
            cells = cells(:, disp_mol);
            i_str = num2str(disp_mol); 
        end
        
        if showI %whether or not to display I (slow down simulation)
            %title(ax, sprintf('t=%d, p = %.2f', t, mean(cells)), ...
            %    'FontSize', 20);
            this_p = mean(cells);
            this_I = moranI(cells, a0*dist);
            title(ax, sprintf('$$t=%d, p%s = %.2f, I%s = %.2f$$', t, i_str, this_p, i_str, this_I), ...
                'FontSize', 20);
        else
            this_p = mean(cells);
            title(ax, sprintf('$$t=%d, p%s = %.2f$$', t, i_str, this_p), ...
                'FontSize', 20);
        end
        %}
        c_all = repmat(1-cells, 1, 3);  
        %disp(c_all);
    end
    % --update cells--
    set(plothandle, 'CData', c_all);
        % to add: different cell types

    % pause before plotting next frame
    pauseTime = 1/10;%/app.SpeedSlider.Value;
    pause(pauseTime);
end