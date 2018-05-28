function [] = tune_B_plot_work_vs_deltahB(B_all, pini, pend, W, qsave, protocol, name)
% Usefulness of code questionable
DeltahB = - B_all(end)*(2*pend-1) + B_all(1)*(2*pini-1);
h = figure();
line = linspace(-4,2,100);
hold on
plot(line, line, 'LineWidth', 1.5); 
scatter(W, DeltahB, 100);
%xlim([-6 2]);
%ylim([-6 0]);
set(gca,'FontSize', 20)
xlabel('$$W$$','FontSize', 24);
ylabel('$$\Delta h_B$$', 'FontSize', 24);
legend('W = \Deltah_B', 'Location','nw');
if qsave
    out_file = fullfile(pwd, 'figures', 'work_distribution', ...
        strcat('W_vs_DeltahB_prot',protocol, name)); % filename
    save_figure_pdf(h, 10, 8, out_file);
end

end