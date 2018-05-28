clear variables
close all
clc

% K = 6, Con = 15
fdata = [2.5 1.5 1.5 1 0.5; 2 1.2 0.8 0.5 0.25;...
    1.5 0.8 0.8 0.5 0.2; 1 5/6 0.5 1/3 0.3];
Nvals = [6 11 15 21];
a0vals = [0.5 0.8 1.0 1.2 1.5];
[a0mesh, Nmesh] = meshgrid(a0vals, Nvals);

set(0, 'defaulttextinterpreter', 'latex');
h1 = figure(1);
imagesc(a0vals, Nvals, fdata);
c = colorbar;
set(gca, 'YDir', 'Normal');
xticks(a0vals);
yticks(Nvals);
xlabel('$$a_0$$');
ylabel('$$\sqrt{N}$$');
%%
f = fit([Nmesh(:), a0mesh(:)], fdata(:), 'poly22');
disp(f)
h2 = figure(2);
plot(f,[Nmesh(:), a0mesh(:)], fdata(:));
xlabel('$$N$$');
ylabel('$$a_0$$');
zlabel('$$\phi$$');
set(gca, 'FontSize', 24);
fname = fullfile(pwd, 'rebuttal', 'I_correction', 'phi_vs_N_a0_plot');
save_figure_pdf(h2, 10, 8, fname);
save_figure_eps(h2, 10, 8, fname);

%%
h3 = figure(3);
%plot(a0vals, fdata, '-o');
semilogy(a0vals, fdata, '-o');
xlabel('$$a_0$$');
ylabel('$$\phi_{a_0}$$');

%%
h4 = figure(4);
hold on
plot(Nvals, fdata', '-o');
ylim([0 2.5]);
xlabel('$$\sqrt{N}$$');
ylabel('$$\phi_{N}$$');
p1vals = zeros(numel(a0vals), 1);
p2vals = p1vals;
for i=1:numel(a0vals)    
    f4 = fit(Nvals', fdata(:, i), 'poly1');
    p1vals(i) = f4.p1;
    p2vals(i) = f4.p2;
    plot(Nvals, f4.p1.*Nvals + f4.p2, '--', 'Color', [0.5 0.5 0.5]);
end
%% Save data
close all;
fname_str = 'I_correction_phi_N_a0_K6_Con15';
fname = fullfile(pwd, 'rebuttal',  'I_correction', fname_str);
save(fname)