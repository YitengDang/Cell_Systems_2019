%%
clear all
close all
warning off
%%
% Plots phenotype maps for system with finite Hill coefficient from saved
% data

% Save figures?
qsave = 1;

% Import data
gridsize = 15;
a0 = 1.5;
hill = 2;
initialID = 'uniform'; %ID for initial configuration
%Examples: 'uniform', 'allON', 'singleOFF'
version = 2;

fname_str = strcat(strrep(sprintf('Con_K_map_a0%d_N%d_hill%.1f_hexagonal_',...
     round(10*a0), gridsize^2, hill), '.', 'p'),...
     initialID, '-v', int2str(version));
disp(fname_str);
fname = fullfile(pwd, 'data','phenotype_finite_Hill',...
            strcat(fname_str, '.mat'));
load(fname);
%%
% 
dK = K(2)-K(1);
dCon = Con(2) - Con(1);

% Map 1: <Xi>
set(0, 'defaulttextinterpreter', 'latex');
f=figure(1);
hold on
% Plot cm
colormap(parula);
%colormap(gray);
imagesc(K, Con, xmean);
h=colorbar;
set(gca, 'YDir', 'normal', 'FontSize', 24)
% Plot points
plot(14*ones(1,4), [15 18 21 30], 'rx', 'LineWidth', 2);
plot(14*ones(1,4), [15 18 21 30], 'ro', 'LineWidth', 1);
xlabel('$$K$$');
ylabel('$$C_{ON}$$');
ylabel(h, '$$\langle X_i \rangle$$', 'Interpreter', 'latex');
xlim([1-dK/2 20+dK/2]);
ylim([2-dCon/2 30+dCon/2]);
%%
% Map 2: std(Xi)
f2=figure(2);
colormap(bone);
%colormap(gray);
imagesc(K, Con, xstd);
h=colorbar;
set(gca, 'YDir', 'normal', 'FontSize', 24)
xlabel('$$K$$');
ylabel('$$C_{ON}$$');
ylabel(h, '$$\sigma(X_i)$$', 'Interpreter', 'latex');
%%
% Map 3: <t_final>
f3=figure(3);
colormap(winter);
%colormap(gray);
imagesc(Con, K, tmean);
h=colorbar;
set(gca, 'YDir', 'normal', 'FontSize', 24)
xlabel('$$K$$');
ylabel('$$C_{ON}$$');
ylabel(h, '$$t_{final}$$', 'Interpreter', 'latex');
%% Save maps
if qsave
out_file = fullfile(pwd, 'figures', 'phenotype_finite_Hill', fname_str);
save_figure_pdf(f, 10, 8, strcat(out_file, initialID, '_meanx2'));
%save_figure_pdf(f2, 10, 8, strcat(out_file, initialID, '_stdx'));
%save_figure_pdf(f3, 10, 8, strcat(out_file, initialID, '_meant'));
end