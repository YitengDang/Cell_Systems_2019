function plot_handle = reset_cell_figure(ax, pos, rcell)
    % format plot
    cla(ax);
    
    % calculate figure dimensions
    N = size(pos,1);
    gz = sqrt(N);
    % all sizes in units of pixels
    Sx = 800; %ax.Position(3); %512;
    Sy = (sqrt(3)/2*(gz-1)+2)/(3/2*(gz+1))*Sx;
    a0 = Sx/(sqrt(N)+1); %Lx/(3/2*(gz+1));
    Rcell = rcell*a0;
    set(ax, 'Position', [100 100 Sx Sy]);
    
    % set image properties
    set(gca, 'YTick', [], 'XTick', [], 'Color', [0.8 0.8 0.8]);
    title(gca, 'Simulate dynamics', 'FontSize', 20);
    Lx = 1;
    d = 2*rcell*Lx/(sqrt(N)+1);
    Ly = sqrt(3)/2*Lx;
    xlim([-d Lx+d]);
    ylim([-d Ly+d]);
    
    %% --plot cells--
    hold on
    % colours
    c_all = ones(N, 3);
    clr_k = zeros(N, 3); % black boundaries
    %markers = {'o', 's'};

    plot_handle = scatter(pos(:,1), pos(:,2), Rcell^2, c_all, 'filled', 'o');
    scatter(pos(:,1), pos(:,2), Rcell^2, clr_k, 'o'); % plot cell boundaries
    
    % Plot box outline
    plot([0 Lx], [0 0], 'k--');
    plot([0 Lx], [Ly Ly], 'k--');
    plot([0 0], [0 Ly], 'k--');
    plot([Lx Lx], [0 Ly], 'k--');
    
    % plot round and square cells separately
    %{
    idx0 = (cell_type==0);
    p0 = scatter(ax,pos(idx0,1), pos(idx0,2), Rcell^2, c_all(idx0, :), 'filled', 'o');
    scatter(ax,pos(idx0,1), pos(idx0,2), Rcell^2, clr_k(idx0, :), 'o'); % plot cell boundaries

    idx1 = (cell_type==1);
    p1 = scatter(ax,pos(idx1,1), pos(idx1,2), Rcell^2, c_all(idx1, :), 'filled', 's');
    scatter(ax,pos(idx1,1), pos(idx1,2), Rcell^2, clr_k(idx1, :), 's'); % plot cell boundaries
    %}
    hold off
    %% Add color bar for Xi
    %{
    c = colorbar;
    c.Ticks = 0:0.2:1;
    set(c, 'FontSize', 12);
    %c.Label.String = '$$X_i$$';
    %ylabel(c, '$$X_i$$', 'Interpreter', 'latex', 'FontSize', 16,...
    %    'Rotation', 90);
    caxis([0 1]);
    map = repmat((1:-0.01:0)', 1, 3);
    %map = 'gray';
    colormap(map);
    %}
end