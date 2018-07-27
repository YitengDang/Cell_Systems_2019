clear all 
close all
clc
%%
% Input parameters
Con = [5 5];
Coff = [1 1];
K = [0 3; 3 3];
M_int = [0 1; -1 1];

% conditions
interactions = zeros(2); %0: ON=OFF, 1: ON~=OFF
for i=1:2
    for j=1:2
        if M_int(i,j)==0
            interactions(i,j) = 0;
        else
            interactions(i,j) = 1 + (Con(j) > K(i,j));
        end
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
%% Statistics
X_out_uniq = unique(X_out_ind); % unique out states
n_uniq = numel(X_out_uniq); % # unique states
 % cycles

%% Draw gene circuit
h1 = figure(1);
hold on
G = digraph(abs(M_int));
s = [0 1];
t = [0 0];

edgecolors = M_int(M_int~=0);

plot(G, 'XData', s, 'YData', t, 'ArrowSize', 20, 'EdgeAlpha', 1, ...
    'EdgeCData', edgecolors, 'LineWidth', 3, 'NodeLabel', {},...
    'Marker', 'o', 'MarkerSize', 100, 'NodeColor', 'k');

nLabels = {'1', '2'};
text(s-0.03,t+0.005, nLabels, 'Color', 'w', 'FontSize', 40); % node labels
%text([-0.5 0.5], [0 0], {'a', 'b'}, 'Color', 'g', 'FontSize', 40); % node labels
%
labels = {'11','12'; '21','22'};
eLabels = {}; %{'a','b','c','d'};
for i=1:2
    for j=1:2
        if M_int(i,j)==0
            continue
        end
        if interactions(i,j)==0
            eLabels{end+1} = labels{i,j};
        elseif interactions(i,j)==1
            eLabels{end+1} = labels{i,j};
        end
    end
end
% positions need to be tuned for each interaction
s2 = [-0.5 0.5; 0.5 1.5];
t2 = [0 -0.05; 0.05 0];
%text(s2(M_int~=0), t2(M_int~=0), eLabels, 'Color', 'g', 'FontSize', 24); 
%}

ax = gca;
%axis([-0.6 1.6 -0.5 0.5]);
%colormap(hsv(2));
map = [1, 0, 0 
    0, 0, 1];
colormap(map);
%ax.Visible = 'off';
h1.Color = [1 1 1];

%set(ax, 'Units', 'Inches', 'Position', [0 0 8 4]);
set(h1, 'Units', 'Inches', 'Position', [1 1 8 4]);

%save_figure_pdf(h1, 9, 4, fullfile(save_path, strcat(fname_str, '_circuit')));
%% Draw state diagram as directed graph
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

save_figure_pdf(h2, 9, 8, fullfile(save_path, strcat(fname_str, '_state_diagram')));
%% Own drawing program (not working)
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