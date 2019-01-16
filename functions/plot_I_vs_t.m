function [msg, I_t, Theta_t] = plot_I_vs_t(cells_hist, t0, a0, dist, option, fig_pos)
% option: 1 = plot I, 2 = plot Theta

str_options = {'I(t)', 'Theta(t)'};

if isempty(cells_hist)
    msg = sprintf(' Unable to plot %s; ', str_options{option});
    I_t = [];
    Theta_t = [];
    return
end

N = size(cells_hist{1}, 1);
s = size(cells_hist{1}, 2);

tmax = numel(cells_hist)-1;
I_t = zeros(tmax+1, s);
Theta_t = zeros(tmax+1, s);

for i=1:tmax+1
    for j=1:s
        %cells = cells_hist{i}{j};
        cells = cells_hist{i};
        %[I(i,j), Theta(i,j)] = moranI(cells, a0*dist);
        [I_t(i,j), Theta_t(i,j)] = moranI(cells(:, j), a0*dist);
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
    lw = 0.5; ps = '-';
end
ps = '-';

for i=1:s
    switch option
        case 1
            y = I_t(:,i);
        case 2
            y = Theta_t(:,i);
    end
    clr = plot_clrs(i,:);
    plot(t0:t0+tmax, y, ps, 'LineWidth', lw, 'Color', clr);
end

%set(0, 'DefaultTextInterpreter', 'latex');
xlabel('$$t$$', 'Interpreter', 'latex');
plot_labels = {'$$I(t)$$', '$$\Theta(t)/f_N$$'};
ylabel(plot_labels{option}, 'Interpreter', 'latex');
legend(num2cell(string(1:s)));
set(gca, 'FontSize', 24);
xlim([t0 t0+tmax])
ylim([-0.1 1]);

set(h2, 'Units', 'inches', 'Position', fig_pos);

msg = sprintf('Successfully plotted %s; ', str_options{option});  
end