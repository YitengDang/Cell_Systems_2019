function msg = plot_I_vs_t_moving_cells(cells_hist, t0, a0, pos_hist,...
    option, fig_pos)
% pos_hist: histogram of cell positions
% t0: starting time
% option: 1 = plot I, 2 = plot Theta

str_options = {'I(t)', 'Theta(t)'};

if isempty(cells_hist)
    msg = sprintf(' Unable to plot %s; ', str_options{option});
    return
end

% Input parameters
N = size(cells_hist{1}, 1);
s = size(cells_hist{1}, 2);

% constants for calculating dist
gz = sqrt(N);
Lx = 1;
delx = Lx/gz;
dely = sqrt(3)/2*delx;
Ly = dely*gz;
    
% variables to store
tmax = numel(cells_hist)-1;
I = zeros(tmax+1, s);
Theta = zeros(tmax+1, s);

for i=1:tmax+1
    pos = pos_hist{i};
    dist = calc_dist_periodic(pos(:,1), pos(:,2), Lx, Ly);
    for j=1:s
        %cells = cells_hist{i}{j};
        cells = cells_hist{i};
        %[I(i,j), Theta(i,j)] = moranI(cells, a0*dist);
        [I(i,j), Theta(i,j)] = moranI(cells(:, j), a0*dist);
    end
end
%%        
h2 = figure(1+option);
cla(h2, 'reset');
hold on
plot_clrs = [1 0 0;
             0 0 1];
if length(cells_hist) < 100
    lw = 1.5; ps = 'o-';
elseif length(cells_hist) < 500
    lw = 1; ps = '.-';
else
    lw = 0.5; ps = '.-';
end
          
for i=1:s
    switch option
        case 1
            y = I(:,i);
        case 2
            y = Theta(:,i);
    end
    clr = plot_clrs(i,:);
    plot(t0:t0+tmax, y, ps, 'LineWidth', lw, 'Color', clr);
end

%set(0, 'DefaultTextInterpreter', 'latex');
xlabel('$$t$$', 'Interpreter', 'latex');
plot_labels = {'$$I(t)$$', '$$\frac{\Theta(t)}{\sum_{i,j}f(r_{ij})}$$'};
ylabel(plot_labels{option}, 'Interpreter', 'latex');
legend(num2cell(string(1:s)));
set(gca, 'FontSize', 24);
xlim([t0 t0+tmax])
ylim([-1 1]);

set(h2, 'Units', 'inches', 'Position', fig_pos);

msg = sprintf('Successfully plotted %s; ', str_options{option});  
end