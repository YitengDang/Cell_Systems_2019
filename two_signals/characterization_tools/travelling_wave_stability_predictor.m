% Predicts whether a travelling wave can propagate according to nearest
% neighbour interactions
clear all
close all
set(0,'defaulttextinterpreter', 'latex')
%% Parameters
gz = 15;
N = gz^2;
a0 = 1.5;
rcell = 0.2;
Rcell = rcell*a0;
lambda = [1 1.2];

Con = [18 16];
K = [0 9; 11 1];

% get pos, dist
mcsteps = 0;
[pos, dist] = initial_cells_random_markov_periodic(gz, mcsteps, rcell);

% calculate fN
idx = gz + round(gz/2); % pick cell not at corner of grid
dist_vec = a0*dist(idx,:);
r = dist_vec(dist_vec>0);
fN = zeros(2,1);
fN(1) = sum(sinh(Rcell)*exp((Rcell-r)./lambda(1)).*(lambda(1)./r));
fN(2) = sum(sinh(Rcell)*exp((Rcell-r)./lambda(2)).*(lambda(2)./r));

% calculate fnn
fnn = zeros(1, 2);
rnn = a0;
fnn(1) = sinh(Rcell)*exp((Rcell-rnn)./lambda(1)).*(lambda(1)./rnn) ; % calculate signaling strength
fnn(2) = sinh(Rcell)*exp((Rcell-rnn)./lambda(2)).*(lambda(2)./rnn) ; % calculate signaling strength

fnnn = zeros(2,2);
rnnn = [sqrt(3) 2].*a0;
fnnn(1,:) = sinh(Rcell)*exp((Rcell-rnnn)./lambda(1)).*(lambda(1)./rnnn);
fnnn(2,:) = sinh(Rcell)*exp((Rcell-rnnn)./lambda(2)).*(lambda(2)./rnnn);

% calculate mean-field contribution
z = 6; % coordination number
p = [2/gz 2/gz]; %single band, horizontal / vertical 
%p = [4/gz 4/gz]; %double band, horizontal / double vertical  / diagonal

Y_mf = zeros(1, 2);
Y_mf(1) = (fN(1) - z*fnn(1))*( Con(1)*p(1) + (1-p(1)) );
Y_mf(2) = (fN(2) - z*fnn(2))*( Con(2)*p(2) + (1-p(2)) );

% save folder 
save_folder = 'H:\My Documents\Multicellular automaton\figures\trav_wave_stability';
fname_str_default = strrep(sprintf('N%d_a0_%.1f_rcell_%.1f_lambda12_%.1f_Con_%d_%d_K_%d_%d_%d_%d',...
    N, a0, rcell, lambda(2), Con(1), Con(2), K(1,1), K(1,2), K(2,1), K(2,2)), '.', 'p');
%% Check for a set of parameters
K_12_all = 1:30;
K_21_all = 1:30;
K_22_all = 1:10;
%check_results = zeros(5, numel(K_12_all), numel(K_21_all), numel(K_22_all));
trav_wave_cond_met = zeros(5, numel(K_12_all), numel(K_21_all), numel(K_22_all));

for k=1:numel(K_22_all)
    for i=1:numel(K_12_all)
        for j=1:numel(K_21_all)
            thisK = K;
            thisK(1,2) = K_12_all(i);
            thisK(2,1) = K_21_all(j);
            thisK(2,2) = K_22_all(k);

            % new method
            Xself = [0 0; 0 1; 1 1; 1 0; 0 0]; % fixed
            
            % plane wave
            Xnei = [0 2; 2 4; 4 4; 4 2; 2 0]; 
            cond_set1 = calc_conditions(fnn, Con, thisK, Y_mf, z, Xself, Xnei);
            trav_wave_cond_met(1, i, j, k) = all(cond_set1);
        
            % inward bend
            Xnei = [0 3; 3 5; 5 3; 3 1; 1 0]; 
            cond_set2 = calc_conditions(fnn, Con, thisK, Y_mf, z, Xself, Xnei);
            trav_wave_cond_met(2, i, j, k) = all(cond_set2);
            
            % outward bend
            Xnei = [0 1; 1 3; 3 5; 5 3; 3 0]; 
            cond_set3 = calc_conditions(fnn, Con, thisK, Y_mf, z, Xself, Xnei);
            trav_wave_cond_met(3, i, j, k) = all(cond_set3);

            % old method
            %{
            cond_set1 = [...
                1+2*Con(2)*fnn(2)+4*fnn(2)+Y_mf(2)<thisK(1,2) ...
                Con(2)+4*Con(2)*fnn(2)+2*fnn(2)+Y_mf(2)>thisK(1,2)...
                1+2*Con(1)*fnn(1)+4*fnn(1)+Y_mf(2)<thisK(2,1)...
                Con(1)+4*Con(1)*fnn(1)+2*fnn(1)+Y_mf(2)>thisK(2,1)...
                1+2*Con(2)*fnn(2)+4*fnn(2)+Y_mf(2)>thisK(2,2)...
                1+6*fnn(2)+Y_mf(2)<thisK(2,2)];
            check_results(i, j, 1) = all(cond_set1);


            cond_set2 = [...
                1+3*Con(2)*fnn(2)+3*fnn(2)+Y_mf(2)<thisK(1,2) ...
                Con(2)+3*Con(2)*fnn(2)+3*fnn(2)+Y_mf(2)>thisK(1,2)...
                1+3*Con(1)*fnn(1)+3*fnn(1)+Y_mf(1)<thisK(2,1)...
                Con(1)+5*Con(1)*fnn(1)+fnn(1)+Y_mf(1)>thisK(2,1)...
                1+3*Con(2)*fnn(2)+3*fnn(2)+Y_mf(2)>thisK(2,2)...
                Con(1)+3*Con(1)*fnn(1)+Y_mf(1)>thisK(2,1) || 1+Con(2)*fnn(2)+5*fnn(2)+Y_mf(2)<thisK(2,2)...
                1+Con(1)*fnn(1)+5*fnn(1)+Y_mf(1)>thisK(2,1) || 1+6*fnn(2)+Y_mf(2)<thisK(2,2)];
            check_results(i, j, 2) = all(cond_set2);

            cond_set3 = [...
            1+3*Con(2)*fnn(2)+3*fnn(2)+Y_mf(2)<thisK(1,2) ...
            Con(2)+3*Con(2)*fnn(2)+3*fnn(2)+Y_mf(2)>thisK(1,2)...
            1+Con(1)*fnn(1)+5*fnn(1)+Y_mf(1)<thisK(2,1)...
            Con(1)+3*Con(1)*fnn(1)+3*fnn(1)+Y_mf(1)>thisK(2,1)...
            1+Con(2)*fnn(2)+5*fnn(2)+Y_mf(2)>thisK(2,2)...
            1+3*Con(1)*fnn(1)+3*fnn(1)+Y_mf(2)>thisK(2,1) || 1+6*fnn(2)+Y_mf(2)<thisK(2,2)];
            %disp(cond_set3);
            check_results(i, j, 3) = all(cond_set3);

            % plane waves possible of type 1, 2 and 3?
            cond_set123 = [...
                1+3*Con(2)*fnn(2)+3*fnn(2)<thisK(1,2) ...
                Con(2)+3*Con(2)*fnn(2)+3*fnn(2)>thisK(1,2)...
                1+3*Con(1)*fnn(1)+3*fnn(1)<thisK(2,1)...
                Con(1)+3*Con(1)*fnn(1)+3*fnn(1)>thisK(2,1)...
                1+6*fnn(2)<thisK(2,2)...
                1+Con(2)*fnn(2)+5*fnn(2)>thisK(2,2)];
            check_results(i, j, 4) = all(cond_set123);

            % plane waves possible of type 1 and 2?
            cond_set12 = [...
                1+3*Con(2)*fnn(2)+3*fnn(2)<thisK(1,2) ...
                Con(2)+3*Con(2)*fnn(2)+3*fnn(2)>thisK(1,2)...
                1+3*Con(1)*fnn(1)+3*fnn(1)<thisK(2,1)...
                Con(1)+4*Con(1)*fnn(1)+2*fnn(1)>thisK(2,1)...
                1+6*fnn(2)<thisK(2,2)...
                1+2*Con(2)*fnn(2)+4*fnn(2)>thisK(2,2)...
                Con(1)+3*Con(1)*fnn(1)+3*fnn(1)>thisK(2,1) || ...
                    1+Con(2)*fnn(2)+5*fnn(2)<thisK(2,2)];
            disp(cond_set12);
            check_results(i, j, 5) = all(cond_set12);
            %}
        end        
    end
end
    
disp('Parameter bounds estimates:');
fprintf('%.2f < K12 < %.2f \n', 1+2*Con(2)*fnn(2)+4*fnn(2),...
    Con(2)+4*Con(2)*fnn(2)+2*fnn(2));
fprintf('%.2f < K21 < %.2f \n', 1+2*Con(1)*fnn(1)+4*fnn(1),...
    Con(1)+4*Con(1)*fnn(1)+2*fnn(1));
fprintf('%.2f < K22 < %.2f \n', 1+6*fnn(2),...
    1+2*Con(2)*fnn(2)+4*fnn(2));

%% Save data
%
save_folder = 'H:\My Documents\Multicellular automaton\figures\trav_wave_stability\data';
fname_str_data = sprintf('Trav_wave_stability_analytical_%s_K12_K21_K22_range', fname_str_default);  
fname = fullfile(save_folder, fname_str_data);
save(fname, 'K_12_all', 'K_21_all', 'K_22_all', 'trav_wave_cond_met');
%}
%Trav_wave_stability_nn_N225_a0_1p5_rcell_0p2_lambda12_1p2_Con_18_16_K_0_9_11_1_K12_K21_K22_range
%% Load data
load_folder = 'H:\My Documents\Multicellular automaton\figures\trav_wave_stability\data';
fname_str_data = sprintf('Trav_wave_stability_analytical_%s_K12_K21_K22_range', fname_str_default);  
fname = fullfile(load_folder, fname_str_data);
load(fname, 'K_12_all', 'K_21_all', 'K_22_all', 'trav_wave_cond_met');

%% Plot results (new method)
wave_types_str = {'plane', 'inward bend', 'outward bend'};

% Plot results 1D
h = figure;
hold on
%bar(K_12_all, check_results);
plot(K_12_all, trav_wave_cond_met(1,:,1,1), 'o', 'MarkerSize', 15);
plot(K_12_all, trav_wave_cond_met(2,:,1,1), 'x', 'MarkerSize', 15);
plot(K_12_all, trav_wave_cond_met(3,:,1,1), '+', 'MarkerSize', 15);
%plot(K_12_all, check_results(:,4), 'o', 'MarkerSize', 25);
%plot(K_12_all, check_results(:,5), 'o', 'MarkerSize', 25);

legend(wave_types_str, 'Location', 'ne');
xlabel('$K^{(12)}$');
ylabel('Travelling wave possible?');
set(gca, 'FontSize', 20, 'YTick', [0 1], 'YTickLabels', {'no', 'yes'});
set(h, 'Units', 'Inches', 'Position', [1 1 10 8]);
ylim([-0.2 1.2]);

qsave = 0;
fname_str = sprintf('Trav_wave_stab_nn_MF_%s_vs_K12', fname_str_default);    
fname = fullfile(save_folder, fname_str);
save_figure(h, 10, 8, fname, '.pdf', qsave);

%% Plot results as heatmap in 2D
K22_idx = 2;

%for K22_idx=1:10
    fname_str_default = strrep(sprintf('N%d_a0_%.1f_rcell_%.1f_lambda12_%.1f_Con_%d_%d_K_%d_%d_%d_%d',...
        N, a0, rcell, lambda(2), Con(1), Con(2), K(1,1), K(1,2), K(2,1), K22_idx), '.', 'p');
    h = figure;
    hold on
    [K_21_mesh, K_12_mesh] = meshgrid(K_21_all, K_12_all);
    %temp1 = squeeze(trav_wave_cond_met(1, :,:,K22_idx));
    idx1 = (trav_wave_cond_met(1, :,:,K22_idx)==1);
    scatter(K_12_mesh(idx1), K_21_mesh(idx1), 100, 'o');
    idx2 = (trav_wave_cond_met(2, :,:,K22_idx)==1);
    scatter(K_12_mesh(idx2), K_21_mesh(idx2), 100, 'x');
    idx3 = (trav_wave_cond_met(3, :,:,K22_idx)==1);
    scatter(K_12_mesh(idx3), K_21_mesh(idx3), 100, '+');
    legend(wave_types_str, 'Location', 'ne');
    xlabel('$K^{(12)}$');
    ylabel('$K^{(21)}$');
    title(sprintf('$K^{(22)} = %.1f $', K_22_all(K22_idx) ));
    set(gca, 'FontSize', 20);
    set(h, 'Units', 'Inches', 'Position', [1 1 10 8]);

    xlim([K_12_all(1)-1 K_12_all(end)+1]);
    ylim([K_21_all(1)-1 K_21_all(end)+1]);

    qsave = 0;
    fname_str = sprintf('Trav_wave_stab_nn_MF_%s_vs_K12_K21', fname_str_default);    
    fname = fullfile(save_folder, fname_str);
    save_figure(h, 10, 8, fname, '.pdf', qsave);
    %close all
%end

%% 3D scatter plot
[K_21_mesh, K_12_mesh, K_22_mesh] = meshgrid(K_21_all, K_12_all, K_22_all);

h = figure;
hold on
grid on
type = 3;
colors = {'b', 'r', 'g'};
idx1 = squeeze(trav_wave_cond_met(type, :,:,:)==1);
s1 = scatter3(K_12_mesh(idx1), K_21_mesh(idx1), K_22_mesh(idx1),...
    150, 'o', 'filled', 'MarkerFaceColor', colors{type});
s1.MarkerFaceAlpha = 0.5;
%{
eps = 0.0;
idx1 = squeeze(trav_wave_cond_met(1, :,:,:)==1);
offset = eps*rand(sum(idx1(:)), 1);
s1 = scatter3(K_12_mesh(idx1)+offset, K_21_mesh(idx1)+offset,...
    K_22_mesh(idx1)+offset, 100, 'o', 'filled');
s1.MarkerFaceAlpha = 0.5;
idx2 = (trav_wave_cond_met(2, :,:,:)==1);
offset = eps*rand(sum(idx2(:)), 1);
s2=scatter3(K_12_mesh(idx2)+offset, K_21_mesh(idx2)+offset, K_22_mesh(idx2)+offset, 100, 'd', 'filled');
s2.MarkerFaceAlpha = 0.5;
idx3 = (trav_wave_cond_met(3, :,:,:)==1);
s3=scatter3(K_12_mesh(idx3), K_21_mesh(idx3), K_22_mesh(idx3),100, 'v', 'filled');
s3.MarkerFaceAlpha = 0.5;
legend(wave_types_str);
%}
xlim([K_12_all(1)-1 K_12_all(end)+1]);
ylim([K_21_all(1)-1 K_21_all(end)+1]);
zlim([K_22_all(1)-1 K_22_all(end)+1]);
az = -45; el = 20;
view(az, el);
xlabel('$$K^{(12)}$$');
ylabel('$$K^{(21)}$$');
zlabel('$$K^{(22)}$$');
set(gca, 'FontSize', 20);
set(h, 'Units', 'inches', 'position', [1 1 10 6]);
daspect([1 1 1])

qsave = 1;
fname_str = strrep(sprintf('Trav_wave_stab_nn_MF_%s_3D_vs_K12_K21_K22_scatter_%s',...
    fname_str_default, wave_types_str{type}), ' ', '_');    
fname = fullfile(save_folder, fname_str);
save_figure(h, 10, 6, fname, '.pdf', qsave);
%}

%% functions
function cond_set = calc_conditions(fnn, Con, K, Y_mf, z, Xself, Xnei)
    % calculate the specific conditions for a given input set
    fnn_mat = repmat(fnn, 5, 1);
    Con_mat = repmat(Con, 5, 1);

    Y_self = (Con_mat-1).*Xself + 1;
    Y_nei = fnn_mat.*(Con_mat.*Xnei + (z-Xnei));
    Y_mf_mat = repmat(Y_mf, 5, 1);
    Y_all_mat = Y_self + Y_nei + Y_mf_mat;

    sgn1 = [-1; 1; 1; -1; -1];
    cond1 = (Y_all_mat(:,2) - K(1,2)).*sgn1 > 0;

    sgn2 = [-1 1; -1 1; 1 -1; 1 -1; 1 -1];
    %sgn2 = [-1 -1 1 1 1; ...
    %        1 1 -1 -1 -1];
    cond2a = (Y_all_mat(:,1) - K(2,1)).*sgn2(:, 1) > 0;
    cond2b = (Y_all_mat(:,2) - K(2,2)).*sgn2(:, 2) > 0;    
    cond2 = zeros(5, 1);
    cond2(1:2) = cond2a(1:2) & cond2b(1:2); % AND for first two conditions
    cond2(3:5) = cond2a(3:5) | cond2b(3:5); % OR for last two conditions

    cond_set = cond1 & cond2;
end

%% Plot results (old method)
%{
h = figure;
hold on
%bar(K_12_all, check_results);
plot(K_12_all, check_results(:,1), 'o', 'MarkerSize', 15);
plot(K_12_all, check_results(:,2), 'x', 'MarkerSize', 15);
plot(K_12_all, check_results(:,3), '+', 'MarkerSize', 15);
%plot(K_12_all, check_results(:,4), 'o', 'MarkerSize', 25);
%plot(K_12_all, check_results(:,5), 'o', 'MarkerSize', 25);

legend({'plane', 'inw. bend', 'outw. bend'}, 'Location', 'ne');
xlabel('$K^{(12)}$');
ylabel('Travelling wave possible?');
set(gca, 'FontSize', 20, 'YTick', [0 1], 'YTickLabels', {'no', 'yes'});
set(h, 'Units', 'Inches', 'Position', [1 1 10 8]);
ylim([-0.2 1.2]);
%}
%% Simulation data
%{
K_12_all = 4:25;
check_results = ones(numel(K_12_all), 3);
check_results([1:3 21:22], 1) = 0;
check_results([1:5 20:22], 2) = 0;
check_results([1:5 20:22], 3) = 0;
%}
%% Plane wave conditions (explicit example) 
%{
Xself = [0 0; 0 1; 1 1; 1 0; 0 0]; %[0 0 1 1 0; 0 1 1 0 0];
Xnei = [0 2; 2 4; 4 4; 4 2; 2 0]; %[0 2 4 4 2; 2 4 4 2 0];

fnn_mat = repmat(fnn, 5, 1);
Con_mat = repmat(Con, 5, 1);

Y_self = (Con_mat-1).*Xself + 1;
Y_nei = fnn_mat.*(Con_mat.*Xnei + (z-Xnei));
Y_mf_mat = repmat(Y_mf, 5, 1);
Y_all_mat = Y_self + Y_nei + Y_mf_mat;

sgn1 = [-1; 1; 1; -1; -1];
cond1 = (Y_all_mat(:,2) - K(1,2)).*sgn1 > 0;

sgn2 = [-1 1; -1 1; 1 -1; 1 -1; 1 -1];
%sgn2 = [-1 -1 1 1 1; ...
%        1 1 -1 -1 -1];
cond2a = (Y_all_mat(:,1) - K(2,1)).*sgn2(:, 1) > 0;
cond2b = (Y_all_mat(:,2) - K(2,2)).*sgn2(:, 2) > 0;    
cond2 = zeros(5, 1);
cond2(1:2) = cond2a(1:2) & cond2b(1:2); % AND for first two conditions
cond2(3:5) = cond2a(3:5) | cond2b(3:5); % OR for last two conditions

cond_set1 = cond1 & cond2;
disp(cond_set1);

% original description
cond_set1 = [...
    1+2*Con(2)*fnn(2)+4*fnn(2)+Y_mf(2)<K(1,2) ...
    Con(2)+4*Con(2)*fnn(2)+2*fnn(2)+Y_mf(2)>K(1,2)...
    1+2*Con(1)*fnn(1)+4*fnn(1)+Y_mf(2)<K(2,1)...
    Con(1)+4*Con(1)*fnn(1)+2*fnn(1)+Y_mf(2)>K(2,1)...
    1+2*Con(2)*fnn(2)+4*fnn(2)+Y_mf(2)>K(2,2)...
    1+6*fnn(2)+Y_mf(2)<K(2,2)];

disp(cond_set1)
%% Inward kink conditions
cond_set2 = [...
    1+3*Con(2)*fnn(2)+3*fnn(2)+Y_mf(2)<K(1,2) ...
    Con(2)+3*Con(2)*fnn(2)+3*fnn(2)+Y_mf(2)>K(1,2)...
    1+3*Con(1)*fnn(1)+3*fnn(1)+Y_mf(1)<K(2,1)...
    Con(1)+5*Con(1)*fnn(1)+fnn(1)+Y_mf(1)>K(2,1)...
    1+3*Con(2)*fnn(2)+3*fnn(2)+Y_mf(2)>K(2,2)...
    Con(1)+3*Con(1)*fnn(1)+Y_mf(1)>K(2,1) || 1+Con(2)*fnn(2)+5*fnn(2)+Y_mf(2)<K(2,2)...
    1+Con(1)*fnn(1)+5*fnn(1)+Y_mf(1)>K(2,1) || 1+6*fnn(2)+Y_mf(2)<K(2,2)];

disp(cond_set2)
%}
%% Outward kink conditions
%{
K = [0 10; 4 2];
cond_set3 = [...
    1+3*Con(2)*fnn(2)+3*fnn(2)+Y_mf(2)<K(1,2) ...
    Con(2)+3*Con(2)*fnn(2)+3*fnn(2)+Y_mf(2)>K(1,2)...
    1+Con(1)*fnn(1)+5*fnn(1)+Y_mf(1)<K(2,1)...
    Con(1)+3*Con(1)*fnn(1)+3*fnn(1)+Y_mf(1)>K(2,1)...
    1+Con(2)*fnn(2)+5*fnn(2)+Y_mf(2)>K(2,2)...
    1+3*Con(1)*fnn(1)+3*fnn(1)+Y_mf(2)>K(2,1) || 1+6*fnn(2)+Y_mf(2)<K(2,2)];

disp(cond_set3)
%}