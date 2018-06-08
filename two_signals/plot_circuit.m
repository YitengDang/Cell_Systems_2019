function plot_circuit(M_int, Con, K)
% plots genetic circuit

% conditions
interactions = zeros(2); %0: ON=OFF, 1: ON~=OFF
for i=1:2
    for j=1:2
        interactions(i,j) = (Con(j) > K(i,j));
    end
end

%% Draw gene circuit
h1 = figure(1);
hold on
G = digraph(abs(M_int));

edgecolors = M_int(M_int~=0);

% Automatic layout
plot(G, 'Layout', 'circle', 'ArrowSize', 20, 'EdgeAlpha', 1, ...
    'EdgeCData', edgecolors, 'LineWidth', 3, 'NodeLabel', {},...
    'Marker', 'o', 'MarkerSize', 100, 'NodeColor', 'k');
nLabels = {'1', '2'};
text([-1 1]-0.1, [0 0]+0.01, nLabels, 'Color', 'w', 'FontSize', 40); % node labels

%
labels = {'x','o'};
eLabels = {}; 
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

%}
ax = gca;
%axis([-0.6 1.6 -0.2 0.2]);
%colormap(hsv(2));
map = [1, 0, 0 
    0, 0, 1];
colormap(map);
ax.Visible = 'off';
h1.Color = [1 1 1];

set(ax, 'Units', 'Inches', 'Position', [0 0 9 4]);
set(h1, 'Units', 'Inches', 'Position', [1 1 8 4]);

end