function update_figure_periodic_scatter_cell_growth(h_cells, h_borders, cells, t, disp_mol, showI, a0, dist, rcell_all, a0_px)
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
        clrs = 1-cells;
        c_all = zeros(size(cells, 1), 3); 
        c_all(:, 3) = clrs(:, 1); % signal 1 present -> Turn on blue channel
        % c_all(:, 2) = clrs(:, 2); % signal 2 present -> Turn on green channel
        c_all(:, 1) = clrs(:, 2); % signal 2 present -> Turn on red channel
        idx_both =  clrs(:, 1) &  clrs(:, 2);
        c_all(:, 2) = idx_both; % signals 1 & 2 present -> Turn on green channel
        
        set(h_cells, 'cdata', c_all);
        
        % --update cell sizes--
        h_gca = gca;
        set(h_gca, 'Units', 'points');
        %{
        Lx = 1;
        gz = sqrt( size(cells, 1) );
        d = 2/gz;
        Sx_px = h_gca.Position(3);
        Lx_px = (Lx/(Lx+2*d))*Sx_px;
        a0_px = Lx_px/(gz); %Lx/(3/2*(gz+1));
        %}
        Rcell_px_all = rcell_all.*a0_px;
        set(h_cells, 'sizedata', (2*Rcell_px_all).^2);
        set(h_borders, 'sizedata', (2*Rcell_px_all).^2);
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
        
        % --update cell sizes--
        h_gca = gca;
        set(h_gca, 'Units', 'points');
        %{
        Lx = 1;
        gz = sqrt( size(cells, 1) );
        d = 2/gz;
        Sx_px = h_gca.Position(3);
        Lx_px = (Lx/(Lx+2*d))*Sx_px;
        a0_px = Lx_px/(gz); %Lx/(3/2*(gz+1));
        %}
        Rcell_px_all = rcell_all.*a0_px;
        set(h_cells, 'sizedata', (2*Rcell_px_all).^2);
        set(h_borders, 'sizedata', (2*Rcell_px_all).^2);
    end
end