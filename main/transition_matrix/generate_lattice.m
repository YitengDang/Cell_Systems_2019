function plot_handle = generate_lattice(pos, rcell, hex_lattice)
    if nargin<3
        hex_lattice = 0;
    end
    h = figure;
    ax = gca;
    %cla(ax);
    ax.Box = 'on';
    ax.Color = 0.8*ones(3,1);
    
    % set image properties
    set(ax, 'YTick', [], 'XTick', []);
    hold(ax, 'on')
    
    scale_factor = 1.5; % scale cells of the hexagonal lattice so they don't look too small
     
    N = size(pos,1);
    gz = sqrt(N);
    %% Plot lattice: (1) hexagonal old style, (2) new style box
    if hex_lattice
        % Show hexagonal lattice
        [~, pos] = init_dist_hex(gz, gz);
        xlim(ax, [-1 pos(end,1)+1]);
        ylim(ax, [-1 pos(end,2)+1]);

        %% calculate figure dimensions
        
        Lx = 1024;
        Ly = (sqrt(3)/2*(gz-1)+2)/(3/2*(gz+1))*Lx;
        a0 = Lx/(3/2*(gz+1));
        Rcell = scale_factor*rcell*a0;
        
        set(h, 'Position', [100 100 Lx Ly]);
        
        % --plot cells--
        % colours
        c_all = ones(N, 3);
        clr_k = zeros(N, 3); % black boundaries
        %markers = {'o', 's'};

        plot_handle = scatter(ax, pos(:,1), pos(:,2), Rcell^2, c_all, 'filled', 'o');
        scatter(ax, pos(:,1), pos(:,2), Rcell^2, clr_k, 'o'); % plot cell boundaries
    else 
        d = 2*rcell*1/(gz+1);
        xlim(ax, [-d 1+d]);
        ylim(ax, [-d sqrt(3)/2+d]);

        % calculate figure dimensions
        Sx = 1024;
        Sy = (sqrt(3)/2*(gz-1)+2)/(3/2*(gz+1))*Sx;
        a0 = Sx/(1+2*d)/(gz+1);
        Rcell = rcell*a0;
        
        set(h, 'Position', [100 100 Sx Sy]);
        
        % --plot cells--
        % colours
        c_all = ones(N, 3);
        clr_k = zeros(N, 3); % black boundaries
        %markers = {'o', 's'};

        plot_handle = scatter(ax, pos(:,1), pos(:,2), Rcell^2, c_all, 'filled', 'o');
        scatter(ax, pos(:,1), pos(:,2), Rcell^2, clr_k, 'o'); % plot cell boundaries

        % plot round and square cells separately
        %{
        idx0 = (cell_type==0);
        p0 = scatter(ax,pos(idx0,1), pos(idx0,2), Rcell^2, c_all(idx0, :), 'filled', 'o');
        scatter(ax,pos(idx0,1), pos(idx0,2), Rcell^2, clr_k(idx0, :), 'o'); % plot cell boundaries

        idx1 = (cell_type==1);
        p1 = scatter(ax,pos(idx1,1), pos(idx1,2), Rcell^2, c_all(idx1, :), 'filled', 's');
        scatter(ax,pos(idx1,1), pos(idx1,2), Rcell^2, clr_k(idx1, :), 's'); % plot cell boundaries
        %}

        % Plot box outline
        Lx = 1; Ly = sqrt(3)/2;
        plot(ax,[0 Lx], [0 0], 'k--');
        plot(ax,[0 Lx], [Ly Ly], 'k--');
        plot(ax,[0 0], [0 Ly], 'k--');
        plot(ax,[Lx Lx], [0 Ly], 'k--');
        hold(ax, 'off');
    end
    %% calculate figure dimensions
    %{
    title(ax, 'Simulate dynamics', 'FontSize', 20);
    % set image properties
    set(ax, 'YTick', [], 'XTick', []);
    hold(ax, 'on')
    xlim(ax, [-1 pos(end,1)+1]);
    ylim(ax, [-1 pos(end,2)+1]);

  
    N = size(pos,1);
    gz = sqrt(N);
    Lx = 512;
    Ly = (sqrt(3)/2*(gz-1)+2)/(3/2*(gz+1))*Lx;
    a0 = Lx/(3/2*(gz+1));
    Rcell = 0.4*a0;
    %app.MessagesTextArea.Value = strcat(app.MessagesTextArea.Value, 'Lx = ', num2str(Lx));
    %app.MessagesTextArea.Value = strcat(app.MessagesTextArea.Value, 'a0 = ', num2str(a0));
    
    % -- set figure position
    set(hin, 'Position', [100 100 Lx Ly]);
    
    % --plot cells--
    % colours
    c_all = ones(N, 3);
    clr_k = zeros(N, 3); % black boundaries
    %markers = {'o', 's'};

    plothandle=scatter(ax, pos(:,1), pos(:,2), Rcell^2, c_all, 'filled', 'o');
    scatter(ax, pos(:,1), pos(:,2), Rcell^2, clr_k, 'o'); % plot cell boundaries

    % plot round and square cells separately
    %{
    idx0 = (cell_type==0);
    p0 = scatter(ax,pos(idx0,1), pos(idx0,2), Rcell^2, c_all(idx0, :), 'filled', 'o');
    scatter(ax,pos(idx0,1), pos(idx0,2), Rcell^2, clr_k(idx0, :), 'o'); % plot cell boundaries

    idx1 = (cell_type==1);
    p1 = scatter(ax,pos(idx1,1), pos(idx1,2), Rcell^2, c_all(idx1, :), 'filled', 's');
    scatter(ax,pos(idx1,1), pos(idx1,2), Rcell^2, clr_k(idx1, :), 's'); % plot cell boundaries
    %}

    hold(ax, 'off');
    %{
    % Add color bar for Xi
    c = colorbar(ax);
    c.Ticks = 0:0.2:1;
    set(c, 'FontSize', 12);
    %c.Label.String = '$$X_i$$';
    %ylabel(c, '$$X_i$$', 'Interpreter', 'latex', 'FontSize', 16,...
    %    'Rotation', 90);
    map = repmat((1:-0.01:0)', 1, 3);
    %map = 'gray';
    colormap(ax, map);
    %}
    %}
end