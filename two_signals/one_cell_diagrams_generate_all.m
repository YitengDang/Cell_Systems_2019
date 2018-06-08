%% Loop over conditions to generate and save all diagrams
clear all 
close all
clc
%%
% Input parameters
Con = [30 30];
Coff = [1 1];
K = 20*[1 1; 1 1];

for m1=-1:1
    for m2=-1:1
        for m3=-1:1
            for m4=-1:1
                
                M_int = [m1 m2; m3 m4];
                disp('M_int = ');
                disp(M_int);
                %

                % conditions
                interactions = zeros(2); %0: ON=OFF, 1: ON~=OFF
                for i=1:2
                    for j=1:2
                        interactions(i,j) = (Con(j) > K(i,j));
                    end
                end

                % Filename
                save_path = 'H:\My Documents\Multicellular automaton\figures\two_signals\one_cell_diagrams';
                fname_str = sprintf('M_int_%d_%d_%d_%d_interactions_%d_%d_%d_%d', M_int(1,1), M_int(1,2),...
                    M_int(2,1), M_int(2,2), interactions(1,1), interactions(1,2),...
                    interactions(2,1), interactions(2,2));

                % state transitions
                X_out_sub = cell(2);
                X_out_ind = zeros(4, 1);
                vec_dir_ind = zeros(4, 2);
                A = zeros(4); % graph adjacency matrix
                for X1=0:1
                    for X2=0:1 
                        ind = sub2ind([2 2], X1+1, X2+1);
                        %disp(ind);
                        X = [X1 X2];
                        C = (Con-Coff).*X+Coff;
                        out = (([C; C] - K).*M_int > 0) + (1 - abs(M_int));
                        X_out = prod(out, 2);
                        X_out_2 = sub2ind([2 2], X_out(1)+1, X_out(2)+1);

                        X_out_sub{X1+1,X2+1} = X_out;
                        %X_out_ind(X1+1,X2+1) = sub2ind([2 2], X_out(1)+1, X_out(2)+1);
                        X_out_ind(ind) = X_out_2;

                        % store as MC graph
                        A(ind, X_out_2) = 1;

                        % find vector directions
                        vec_dir_ind(ind, :) = (X_out - X');
                    end
                end
                % Statistics
                %X_out_uniq = unique(X_out_ind); % unique out states
                %n_uniq = numel(X_out_uniq); % # unique states
                % cycles

                %% Draw gene circuit
                plot_circuit(M_int, Con, K);

                %{
                h1 = figure(1);
                hold on
                G = digraph(abs(M_int));

                edgecolors = M_int(M_int~=0);

                % (1) manually setting node positions
                %{
                s = [0 5];
                t = [0 0];

                plot(G, 'XData', s, 'YData', t, 'ArrowSize', 20, 'EdgeAlpha', 1, ...
                    'EdgeCData', edgecolors, 'LineWidth', 3, ...
                    'Marker', 'o', 'MarkerSize', 100, 'NodeColor', 'k');
                nLabels = {'1', '2'};
                text(s-0.03,t+0.005, nLabels, 'Color', 'w', 'FontSize', 40); % node labels
                %}
                %{ 
                % -- only self-interaction
                G1 = digraph([1]);
                G2 = digraph([1]);
                edgecolors = [0 1];
                plot(G1,'XData', s(1), 'YData', t(1), 'ArrowSize', 20, 'EdgeAlpha', 1, ...
                    'EdgeCData', edgecolors(1), 'LineWidth', 3, ...
                    'Marker', 'o', 'MarkerSize', 100, 'NodeColor', 'k');
                plot(G2,'XData', s(2), 'YData', t(2), 'ArrowSize', 20, 'EdgeAlpha', 1, ...
                    'EdgeCData', edgecolors(2), 'LineWidth', 3, ...
                    'Marker', 'o', 'MarkerSize', 100, 'NodeColor', 'k');
                %}
                % (2) Automatic layout
                plot(G, 'Layout', 'circle', 'ArrowSize', 20, 'EdgeAlpha', 1, ...
                    'EdgeCData', edgecolors, 'LineWidth', 3, 'NodeLabel', {},...
                    'Marker', 'o', 'MarkerSize', 100, 'NodeColor', 'k');
                nLabels = {'1', '2'};
                text([-1 1]-0.1, [0 0]+0.01, nLabels, 'Color', 'w', 'FontSize', 40); % node labels

                labels = {'x','o'};
                eLabels = {}; %{'a','b','c','d'};
                for i=1:2
                    for j=1:2
                        if M_int(i,j)==0
                            continue
                        end
                        if interactions(i,j)==0
                            eLabels{end+1} = labels{1};
                        elseif interactions(i,j)==1
                            eLabels{end+1} = labels{2};
                        end
                    end
                end
                s2 = [-2.2 -0.05; -0.05 2.1];
                t2 = [0.01 -0.15; 0.15 0.01];
                text(s2(M_int~=0), t2(M_int~=0), eLabels, 'Color', [0 0.7 0], 'FontSize', 24); 

                ax = gca;
                %axis([-2.3 2.3 -0.6 0.6]);
                %colormap(hsv(2));
                map = [1, 0, 0 
                    0, 0, 1];
                colormap(map);
                ax.Visible = 'off';
                h1.Color = [1 1 1];

                set(ax, 'Units', 'Inches', 'Position', [0 0 9 4]);
                set(h1, 'Units', 'Inches', 'Position', [1 1 8 4]);
                %}
                fname = fullfile(save_path, strcat(fname_str, '_circuit'));
                if exist(strcat(fname, '.pdf'), 'file')~=2
                    save_figure_pdf(gcf, 9, 4, fname);
                end
                %% Draw state diagram as directed graph
                plot_state_diagram(M_int, Con, Coff, K);

                %{
                h2 = figure(2);
                hold on
                s = [0 1 0 1];
                t = [0 0 1 1];
                %A = ones(4);
                Gs = digraph(A);
                nLabels = {};
                plot(Gs, 'XData', s, 'YData', t, 'ArrowSize', 20, 'EdgeAlpha', 1, ...
                    'LineWidth', 3, 'EdgeColor', 'k',...
                    'Marker', 'o', 'MarkerSize', 100, 'NodeColor', [0.2 0.2 0.2], 'NodeLabel', nLabels);
                text(s-0.09,t+0.015,{'(0,0)','(0,1)','(1,0)','(1,1)'}, 'Color','w', 'FontSize', 32)
                ax = gca;
                axis([-0.4 1.4 -0.4 1.4]);
                ax.Visible = 'off';
                h2.Color = [1 1 1];
                set(ax, 'Units', 'Inches', 'Position', [0 0 9 8]);
                set(h2, 'Units', 'Inches', 'Position', [1 1 9 8]);
                %}
                fname = fullfile(save_path, strcat(fname_str, '_state_diagram'));
                if exist(strcat(fname, '.pdf'), 'file')~=2
                    save_figure_pdf(gcf, 9, 8, fname);
                end
                
                %}
                %%
                close all
            end
        end
    end
end
               
%% Manual drawing program (not working)
%{
h=figure(1);
set(gcf, 'Units', 'Inches', 'Position', [1 1 9 8]);

hold on
axis([0 1.2 0 1.2]);
xpos = [0 1];
ypos = [0 1];
xw = 0.2; %width
yw = 0.2;
for i=1:2
    for j=1:2
        ind = sub2ind([2 2], i, j);
        pos = [xpos(i) ypos(j) xw yw];
        rectangle('Position', pos, 'Curvature', [0 0]); 
        
        quiver(xpos(i)+xw/2, ypos(j)+xw/2, vec_dir_ind(ind, 1), ...
            vec_dir_ind(ind, 2), 'b'); % transition arrow
    end
end

drawnow;
%}