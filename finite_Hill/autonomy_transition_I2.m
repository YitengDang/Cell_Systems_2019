% Plot average value of final I2 (modified spatial order parameter) for
% different values of a0
clear all
close all
%warning off
set(0, 'defaulttextinterpreter', 'latex');

%%
% lattice parameters
gridsize = 11;
N = gridsize^2;
%a0 = 5.2;
%Rcell = 0.2*a0;

% circuit parameters
Con = 16;
K = 8;
hill = 2;

% initial conditions
initialID = 'binaryrand'; %'binary'
p0 = 0.3;
iniON = round(p0*N);

[dist, pos] = init_dist_hex(gridsize, gridsize);
%% Loop over a0
a0list = 5.25:0.005:5.30; %5.0:0.05:5.9;
I2mean = zeros(numel(a0list), 1);
I2std = I2mean;
I2mean_ord = I2mean;
I2std_ord = I2mean;
I2ordcount = I2mean;
nonunif_frac = zeros(numel(a0list), 1);

% histogram
edges = -0.08:0.02:0.18;
I2_hist = zeros(numel(a0list), numel(edges)-1);

for idx = 1:numel(a0list)
a0 = a0list(idx);
%Rcell = 0.2*a0;
%dist_vec = a0*dist(1,:);
%r = dist_vec(dist_vec>0); % exclude self influence
%fN = sum(sinh(Rcell)*sum(exp(Rcell-r)./r)); % calculate signaling strength
%% Load data, analyze and plot
% Path to search for the saved data. It searchs by the name, defined by the
% parameters chosen
path = fullfile(pwd, 'data', 'dynamics_transition 2017-10-25'); 
straux = '(\d+)';
straux2 = '(\w+)';
if a0==5.25 || a0==5.30
    fpattern = strrep(sprintf('N%d_a0_%.2f_Con%.2f_K%.2f_hill%.2f_%s-v%s', ...
        N, a0, Con, K, hill, initialID, straux), '.', 'p');
else
    fpattern = strrep(sprintf('N%d_a0_%.3f_Con%.2f_K%.2f_hill%.2f_%s-v%s', ...
        N, a0, Con, K, hill, initialID, straux), '.', 'p');    
end

% Get all file names in the directory
listing = dir(path);
num_files = numel(listing)-2; %first two entries are not useful
count = 0;
for i = 1:num_files
    filename = listing(i+2).name;
    % remove extension and do not include txt files
    [~,name,ext] = fileparts(filename);
    if strcmp(ext, '.mat')
        count = count + 1;
        names{count} = name;
    end
end

% Plot all final states
%h1 = figure(1);
%hold on

% Save I2 final
I2fin = [];
nonunif_count = 0;
count = 0;
for i = 1:numel(names)
    [tokens, ~] = regexp(names{i},fpattern,'tokens','match');
    if numel(tokens) > 0
        disp(names{i}) % displays the file name
        count = count+1;
        % load the data
        load(fullfile(path,strcat(names{i},'.mat')));
        %cells_final = cells_hist{end};
        
        [~, theta] = moranI(cells_final, a0*dist);
        p = sum(cells_final)/N;
        I2fin(end+1) = theta - (2*p-1)^2;
        
        %figure(h1); 
        %plot(repmat(count, N, 1), cells_final, 'x');
        
        % count final states which are nonuniform
        %{
        if sum(round(cells_final))==0 || sum(round(cells_final))==N
            nonunif_count = nonunif_count + 1;
        else
            % show nonuniform lattice
            hin = figure();
            cell_type = zeros(N, 1);
            
            update_cell_figure_continuum(hin, pos, a0, round(cells_final), cell_type, 0)
            %k = waitforbuttonpress;
        end
        %}
    end
end
nonunif_frac(idx) = nonunif_count/count;

% Heat map of I2 values
[I2_hist(idx, :), edges] = histcounts(I2fin, edges);


% Get mean and std of I2 values
I2mean(idx) = mean(I2fin);
I2std(idx) = std(I2fin);

% Ordered states: define as having sigma < eps
eps = 10^(-9);
I2mean_ord(idx) = mean(I2fin(abs(I2fin)>eps));
I2std_ord(idx) = std(I2fin(abs(I2fin)>eps));
I2ordcount(idx) = sum(abs(I2fin)>eps);
end

%% Plot mean I2 against a0
h100=figure(100);
errorbar(a0list, I2mean, I2std, 'LineWidth', 2);
xlabel('$$a_0$$');
ylabel('$$\langle I_2(t_f) \rangle$$');
set(gca,'FontSize', 24);
set(gcf, 'Units', 'Inches', 'Position', [0.5 0.5 10 8]);
%xlim([a0list(1)-0.05 a0list(end)+0.05]);
xlim([a0list(1)-0.005 a0list(end)+0.005]);
%ylim([0 70]);

qsave = 0;
if qsave
    fname_str = strrep(sprintf('I2_final_avg_a0_%.2fto%.2f_N%d_K%d_Con%d_hill_%.1f_%s',...
        a0list(1), a0list(end), N, K, Con, hill, initialID), '.', 'p');
    fname = fullfile(pwd, 'figures', 'autonomy_transition_I2_final', fname_str);
    save_figure_pdf(h100, 10, 8, fname);
    save_figure_eps(h100, 10, 8, fname);
end

%% Number of ordered final states
h102 = figure(102);
plot(a0list, I2ordcount, '-o');
xlabel('a0');
ylabel('Number of ordered final states');
%}

%% Heat map
h103 = figure(103);
bincenters = (edges(1:end-1)+edges(2:end))/2;
nsim = numel(I2fin); % number of simulations
imagesc(a0list, bincenters, I2_hist'/nsim);
c = colorbar;
xlabel('$$a_0$$');
ylabel('$$P(I_2|a_0)$$');
ylabel(c, 'Probability');
set(gca,'YDir', 'Normal', 'FontSize', 24);
caxis([0 1]);

qsave = 0;
if qsave
    fname_str = strrep(sprintf('I2_final_heat_map_a0_%.2fto%.2f_N%d_K%d_Con%d_hill_%.1f_%s',...
        a0list(1), a0list(end), N, K, Con, hill, initialID), '.', 'p');
    fname = fullfile(pwd, 'figures', 'autonomy_transition_I2_final', fname_str);
    save_figure_pdf(h103, 14, 8, fname);
    save_figure_eps(h103, 14, 8, fname);
end