%% Plot and save all one signal diagrams
clear all 
close all
%%
M_int_all = [1 -1];
texts = {'', 'U', 'P1', 'A1', 'A0', 'P0', 'A01'};
for M_idx = 1:1 % act./deact.
    for phase = 1:1
        M_int = M_int_all(M_idx);
        % Get adjacency matrix
        A_map = cell(2, 6);
        % activation 
        A_map{1,1} = ones(2);
        A_map{1,2} = [0 1; 0 1];
        A_map{1,3} = [1 1; 0 1];
        A_map{1,4} = [1 0; 1 1];
        A_map{1,5} = [1 0; 1 0];
        A_map{1,6} = [1 0; 0 1];
        % repression
        A_map{2,1} = ones(2);
        A_map{2,2} = [1 0; 1 0];
        A_map{2,3} = [1 1; 1 0];
        A_map{2,4} = [0 1; 1 1];
        A_map{2,5} = [0 1; 0 1];
        A_map{2,6} = [0 1; 1 0];

        A = A_map{M_idx, phase};
        %% Draw state diagram
        
        t = texts{phase};
        
        h = figure(10);
        plot(1, 1);
        text(0,0, t, 'FontSize', 160, 'HorizontalAlignment', 'center');
        axis([-0.4 0.4 -0.03 0.03]);
        ax = gca;
        ax.Visible = 'off';
        %{
        hold on
        s = [0 1];
        t = [0 0];
        Gs = digraph(A);
        nLabels = {};
        g=plot(Gs, 'XData', s, 'YData', t, 'ArrowSize', 20, 'EdgeAlpha', 1, ...
            'LineWidth', 3, 'EdgeColor', 'k', ...
            'Marker', 'o', 'MarkerSize', 100, 'NodeColor', [0.2 0.2 0.2], 'NodeLabel', nLabels);
        % Make edges dashed if state has two outgoing edges
        for i=1:2
            if sum(A(i,:))==2
                idx = find(A(i,:));
                highlight(g, i, idx, 'LineStyle', '--', 'LineWidth', 2);
            end
        end
        text(s-0.03,t,{'0', '1'}, 'Color', 'w', 'FontSize', 32)
        ax = gca;
        if A(1,2)==0 && A(2,1)==0
            axis([-0.4 1.4 -0.03 0.03]);
        else
            axis([-0.4 1.4 -0.25 0.25]);
        end
        ax.Visible = 'off';
        h.Color = [1 1 1];
        %set(ax, 'Units', 'Inches', 'Position', [0 0 9 8]);
        %set(h2, 'Units', 'Inches', 'Position', [1 1 9 8]);
        set(ax, 'Units', 'Inches', 'Position', [0 0 7 4]);
        set(h, 'Units', 'Inches', 'Position', [0.2 0.2 7 4]);
        %}
        % Save figures
        %save_folder = 'H:\My Documents\Multicellular automaton\latex\12_common_circuits\figures\all_topologies_one_signal';
        save_folder = 'H:\My Documents\Multicellular automaton\figures\two_signals\all_topologies\all_diagrams_by_topology\text_labels';
        %fname_str = sprintf('state_diagram_M_int_%d_phase_%d', M_int, phase);
        fname_str = sprintf('text_nothing', phase);
        fname = fullfile(save_folder, fname_str);
        save_figure(h, 4, 4, fname, '.pdf', 1);
        
        %pause(1);
        close all;
    end
end
