% Plots the equal-time correlation function averaged over trajectories

clear variables
close all
warning off

% Parameters of the system
gridsize = 11;
N = gridsize^2;
a0 = 0.5;
Rcell = 0.2*a0;
p_initial = 0.8;
iniON = round(p_initial*N);
% parameters of the circuit
Con = 8;
K = 16;

% array to store correlation function values in
tmax = 12; % maximum time from the data set
cells = zeros(N,1);
[dist, ~] = init_dist_hex(gridsize, gridsize);
[~, distvals] = correlation_func(dist, cells);
corr_sum = zeros(numel(distvals),tmax);
corr_counts = zeros(numel(distvals),tmax);
%% load data
path = 'H:\My Documents\Multicellular automaton\data\dynamics\no_noise'; 
straux = '(\d+)';
fpattern = sprintf('N%d_n%d_neq_%s_a0%d_K_%d_Son_%d_t_%s-v%s', ...
    N, iniON, straux, 10*a0, K, Con, straux, straux);

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

for i = 1:numel(names)
    % first get the filename given by the student
    [tokens, ~] = regexp(names{i},fpattern,'tokens','match');
    if numel(tokens) > 0
        disp(names{i}) % displays the file name
        % load the data
        load(fullfile(path,strcat(names{i},'.mat')), 'corr');
        t1 = length(corr);
        for t=1:t1
            corr_sum(:, t) = corr_sum(:, t) + corr{t};
            corr_counts(:, t) = corr_counts(:, t)+1;
        end
    end
end

%% Plot correlation function
corr_func = corr_sum./corr_counts;
h=figure();
set(gca, 'FontSize', 20);
xlim([0 max(distvals)]);
ylim([-0.05 0.25]);    
xlabel('Distance');
ylabel('Mean correlation');
hold on
for i=1:tmax
    plot(distvals, corr_func(:,i), 'o-');
    %pause(0.5);
end
fname_str = sprintf('corr_N%d_iniON%d_a0_%d_K%d_Con%d_nruns%d_tmax%d',...
    N, iniON, 10*a0, K, Con, corr_counts(1,1), tmax);
out_file = fullfile(pwd,'figures','correlation_function', fname_str);
save_figure_pdf(h, 10, 8, out_file);

%% Create GIF animations of figures
options = {'single', 'all'}; %plot functions together or separate
for k=1:2
    opt = options{k};
    h=figure();
    if strcmp(opt, 'all')
        hold on
    end
    axis tight manual % this ensures that getframe() returns a consistent size

    fname_str = sprintf('corr_N%d_iniON%d_a0_%d_K%d_Con%d_tmax%d_nruns%d_%s.gif',...
        N, iniON, 10*a0, K, Con, tmax, corr_counts(1,1), opt);
    out_file = fullfile(pwd,'figures','correlation_function', fname_str);
    for t=1:tmax
        plot(distvals, corr_func(:,t), 'o-');
        set(gca, 'FontSize', 16);
        title(sprintf('t=%d', t));
        xlim([0 max(distvals)]);
        ylim([-0.05 0.25]);
        xlabel('Distance');
        ylabel('Mean correlation');
        drawnow

        % Capture the plot as an image 
        frame = getframe(h); 
        im = frame2im(frame); 
        [imind,cm] = rgb2ind(im,256); 

        % Write to the GIF File 
        if t == 1 
          imwrite(imind,cm,out_file,'gif', 'Loopcount',inf); 
        else 
          imwrite(imind,cm,out_file,'gif', 'WriteMode','append'); 
        end 
    end

    % Not working figure size
    %set(gcf, 'PaperPositionMode', 'manual');
    %set(gcf, 'PaperUnits', 'inches');
    %set(gcf, 'PaperSize', [15 11]);
    %set(gca, 'Units', 'inches');
    %set(gca, 'Position', [0 0 10 8]);
end