function [h1, msg, p_t] = plot_p_vs_t(cells_hist, t0, fig_pos)
    if isempty(cells_hist)
        msg = ' Unable to plot p(t); ';
        p_t = [];
        h1 = [];
        return
    end

    s = size(cells_hist{1}, 2);

    N = size(cells_hist{1}, 1);
    tmax = numel(cells_hist)-1;
    Non = zeros(tmax+1, s);

    for i=1:tmax+1
        cells = cells_hist{i};
        for j=1:s
            Non(i,j) = sum(cells(:, j));
        end
    end
    
    h1 = figure(1);
    cla(h1, 'reset');
    hold on
    plot_clrs = [1 0 0;
                0 0 1];
    if length(cells_hist) < 100
        ps = '-'; lw = 1.2;
    elseif length(cells_hist) < 500
        ps = '.-'; lw = 1;
    else
        ps = '-'; lw = 0.5;
    end
    for i=1:s
        clr = plot_clrs(i, :);
        p1 = plot(t0:t0+tmax, Non(:,i)/N, ps, 'LineWidth', 2, 'Color', clr);
        p1.Color(4) = 0.75; % transparency
    end
    xlabel('t', 'Interpreter', 'tex');
    ylabel('p', 'Interpreter', 'tex');
    legend(num2cell(string(1:s)));
    set(gca, 'FontSize', 24);
    xlim([t0 t0+tmax])
    ylim([0 1]);

    set(h1, 'Units', 'inches', 'Position', fig_pos);

    msg = 'Successfully plotted p(t); ';    

    % Output
    p_t = Non/N; %p(t) 

end