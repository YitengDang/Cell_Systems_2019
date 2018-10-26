% Plots the state diagram required for a travelling wave in terms of E, F,
% M and B states
clear all
close all
clc
%%
A = zeros(4);
A(1,1) = 1;
A(1,2) = 1;
A(2,4) = 1;
A(3,1) = 1;
A(4,3) = 1;
h = draw_state_diagram(A, 1, 0);

subdiagrams{1} = [1 2; 2 4; 3 1; 4 3]; % sets of (i,j) indices of transitions that need to be present

qsave = 0;
folder = 'H:\My Documents\Multicellular automaton\latex\13_pattern_analysis\figures';
fname_str = 'travelling_wave_state_diagram_FMBE_states';
fname = fullfile(folder, fname_str);
save_figure(h, 0, 0, fname, '.pdf', qsave);

function h = draw_state_diagram(A, fig_num, count)
    % A: state diagram (graph adjacency matrix)
    % fig_num: number of figure to plot
    % count: number of times this diagram appears among all phases
    % h: returns figure handle
    if nargin<2
        fig_num = 1;
        count = 0;
    end
    h = figure(fig_num);
    hold on
    s = [0 1 0 1];
    t = [0 0 1 1];
    %A = ones(4);
    Gs = digraph(A);
    nLabels = {};
    g=plot(Gs, 'XData', s, 'YData', t, 'ArrowSize', 20, 'EdgeAlpha', 1, ...
        'LineWidth', 3, 'EdgeColor', 'k',...
        'Marker', 'o', 'MarkerSize', 100, 'NodeColor', [0.2 0.2 0.2], 'NodeLabel', nLabels);
    % Make edges dashed if state has two outgoing edges
    for i=1:4
        if sum(A(i,:))==2
            idx = find(A(i,:));
            highlight(g, i, idx, 'LineStyle', '--', 'LineWidth', 2);
        elseif sum(A(i,:))==4
            highlight(g, i, 1:4, 'LineStyle', ':', 'LineWidth', 2);
        end
    end
    %text(s-0.11,t+0.019,{'(0,0)','(1,0)','(0,1)','(1,1)'}, 'Color', 'w', 'FontSize', 32)
    text(s-0.05,t+0.01,{'E','F','B','M'}, 'Color', 'w', 'FontSize', 40)
    
    %text(0.5, 1.3, sprintf('Number of diagrams = %d', count), ...
    %    'FontSize', 16, 'HorizontalAlignment', 'center')
    ax = gca;
    axis([-0.4 1.4 -0.4 1.4]);
    ax.Visible = 'off';
    h.Color = [1 1 1];
    set(ax, 'Units', 'Inches', 'Position', [0 0 7 6]);
    set(h, 'Units', 'Inches', 'Position', [0.2 0.2 7 6]);
end