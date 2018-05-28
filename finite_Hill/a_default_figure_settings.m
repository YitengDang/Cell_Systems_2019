%% basic properties
set(0, 'defaulttextinterpreter', 'latex');
set(gca, 'FontSize', 24);
set(gcf, 'Units', 'Inches', 'Position', [1 1 10 8]);

%% figures with colorbar
colormap('summer');
c = colorbar;
ylabel(c, 'h=H/N', 'FontSize', 24)

%% saving figures
h = figure();
qsave = 0;
if qsave
    fname_out = 'example'
    out_file = fullfile(pwd, 'temp', fname_out); % filename
    save_figure_pdf(h, 10, 8, out_file);
    save_figure_eps(h, 10, 8, out_file);
end

%% Creating a legend from numeric array
legendCell = cellstr(num2str([1:5]', 'N=%-d'));

%% Set log axis on plot
figure();
p1 = plot(0:10, (0:10).^2);
set(get(p1,'Parent'), 'XScale', 'log');