clear all 
close all
clc
%%
% Input parameters
Con = [15 12];
Coff = [1 1];
K = [50 5; 5 2];
M_int = [1 1; 1 -1];

% conditions
interactions = zeros(2,2);
for i=1:2
    for j=1:2
        interactions(i,j) = (Con(j) > K(i,j));
    end
end

% state transitions
X_out_sub = cell(2);
X_out_ind = zeros(4, 1);
%vec_dir_ind = zeros(4, 2);
A = zeros(4); % graph adjacency matrix
for X1=0:1
    for X2=0:1 
        ind = sub2ind([2 2], X1+1, X2+1);
        %disp(ind);
        X = [X1 X2];
        C = (Con-Coff).*X+Coff;
        out = (([C; C] - K).*M_int > 0) + (1 - abs(M_int));
        X_out = prod(out, 2);
        % if no connections to a gene, output = input (remains constant)
        X_out(sum(abs(M_int), 2)==0) = X(sum(abs(M_int), 2)==0);
        
        X_out_2 = sub2ind([2 2], X_out(1)+1, X_out(2)+1);
        
        X_out_sub{X1+1,X2+1} = X_out;
        %X_out_ind(X1+1,X2+1) = sub2ind([2 2], X_out(1)+1, X_out(2)+1);
        X_out_ind(ind) = X_out_2;
        
        % store as MC graph
        A(ind, X_out_2) = 1;
        
        % find vector directions
        %vec_dir_ind(ind, :) = (X_out - X');
    end
end
%%
%{
X = [0 1];
C = (Con-Coff).*X+Coff;
out = (([C; C] - K).*M_int > 0) + (1 - abs(M_int));
X_out = prod(out, 2);
% if no connections to a gene, output = input (remains constant)
X_out(sum(abs(M_int), 2)==0) = X(sum(abs(M_int), 2)==0);

%% Statistics
X_out_uniq = unique(X_out_ind); % unique out states
n_uniq = numel(X_out_uniq); % # unique states
 % cycles
%}
%% Draw gene circuit
h1 = figure(1);
hold on
M2 = M_int'; 
G = digraph(abs(M2));
s = [-1 1];
t = [0 0];

edgecolors = M2(M2~=0);

% (1) automatic layout
%plot(G, 'Layout', 'circle', 'ArrowSize', 20, 'EdgeAlpha', 1, ...
%    'EdgeCData', edgecolors, 'LineWidth', 3, 'NodeLabel', {},...
%    'Marker', 'o', 'MarkerSize', 100, 'NodeColor', 'k');

% (2) manual layout
g = plot(G, 'XData', s, 'YData', t, 'ArrowSize', 20, 'EdgeAlpha', 1, ...
    'EdgeCData', edgecolors, 'LineWidth', 3, ...
    'Marker', 'o', 'MarkerSize', 100, 'NodeColor', 'k', 'NodeLabel', {});


% Labels
%{
text(s-0.03,t+0.005,{'1', '2'}, 'Color', 'w', 'FontSize', 40)
eLabelList = {'x', 'o'};
eLabels = {};
for i=1:2
    for j=1:2
        if M2(i,j)==0
            continue
        end
        if interactions(j,i)==0 %NB transposed!
            eLabels{end+1} = eLabelList{1};
        elseif interactions(j,i)==1
            eLabels{end+1} = eLabelList{2};
        end
    end
end
s2 = [-2.2 -0.05; -0.05 2.1];
t2 = [0.01 -0.15; 0.15 0.01];
text(s2(M2~=0), t2(M2~=0), eLabels, 'Color', [0 0.7 0], 'FontSize', 24); 
%}
% Line styles
% solid lines -> interaction threshold reached
% dashed lines -> interaction threshold not reached
for i=1:2
    for j=1:2
        if ~interactions(i,j) && M_int(i,j)~=0
            highlight(g, i, j, 'LineStyle', ':')
        end
    end
end

% Style
ax = gca;
%axis([-1 3 -0.2 0.2]);
%colormap(hsv(2));
map = [0, 0, 1
    1, 0, 0];
colormap(map);
%ax.Visible = 'off';
h1.Color = [1 1 1];

%set(ax, 'Units', 'Inches', 'Position', [0 0 8 4]);
set(h1, 'Units', 'Inches', 'Position', [1 1 8 4]);

% save
qsave = 0;
if qsave
    save_path = 'H:\My Documents\Multicellular automaton\figures\two_signals\one_cell_diagrams';
    fname_str = sprintf('M_int_%d_%d_%d_%d_interactions_%d_%d_%d_%d', M_int(1,1), M_int(1,2),...
        M_int(2,1), M_int(2,2), interactions(1,1), interactions(1,2),...
        interactions(2,1), interactions(2,2));
    fname = fullfile(save_path, strcat(fname_str, '_circuit'));
    if exist(strcat(fname, '.pdf'), 'file')~=2
        save_figure_pdf(gcf, 9, 4, fname);
    end
end
%% Draw state diagram as directed graph
h2 = figure(2);
hold on
s = [0 1 0 1];
t = [0 0 1 1];
%A = ones(4);
Gs = digraph(A);
nlabels = {};
plot(Gs, 'XData', s, 'YData', t, 'ArrowSize', 20, 'EdgeAlpha', 1, ...
    'LineWidth', 3, 'EdgeColor', 'k',...
    'Marker', 'o', 'MarkerSize', 100, 'NodeColor', [0.2 0.2 0.2], 'NodeLabel', nlabels);
text(s-0.09,t+0.015,{'(0,0)','(1,0)','(0,1)','(1,1)'}, 'Color','w', 'FontSize', 30)
ax = gca;
axis([-0.4 1.4 -0.4 1.4]);
ax.Visible = 'off';
h2.Color = [1 1 1];
set(ax, 'Units', 'Inches', 'Position', [0 0 9 8]);
set(h2, 'Units', 'Inches', 'Position', [1 1 9 8]);

% save
qsave = 0;
if qsave
    save_path = 'H:\My Documents\Multicellular automaton\figures\two_signals\one_cell_diagrams';
    fname_str = sprintf('M_int_%d_%d_%d_%d_interactions_%d_%d_%d_%d', M_int(1,1), M_int(1,2),...
        M_int(2,1), M_int(2,2), interactions(1,1), interactions(1,2),...
        interactions(2,1), interactions(2,2));
    fname = fullfile(save_path, strcat(fname_str, '_state_diagram'));
    if exist(strcat(fname, '.pdf'), 'file')~=2
        save_figure_pdf(gcf, 9, 8, fname);
    end
end
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