% Plot all parameters of a saved dynamics data for all the files saved that
% have the set of parameters established.

% This plots the lines in p vs I space for several saved runs on top of the
% hamiltonian map. It also plots the history of the hamiltonian
clear variables
close all
warning off

% Parameters of the system
gridsize = 11;
N = gridsize^2;
a0list = 5.4;
%Rcell = 0.2*a0;
initialID = 'binaryrand';

% parameters of the circuit
K = 8;
Con = 16;
hill = 2;
prec = 8;
tmax = 3000;

% Path to search for the saved data. It searchs by the name, defined by the
% parameters chosen
%path = 'H:\My Documents\Multicellular automaton\finite_Hill\data\dynamics 2017-10-16'; 
path = 'H:\My Documents\Multicellular automaton\finite_Hill\data\dynamics binaryrand';
straux = '(\d+)';
straux2 = '(\w+)';

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

%{
% Find uniform lattice fixed points
fp = zeros(3, 1);
x0 = [0.03 0.2 0.65]; %estimates based on previous graph
%hfunc = @update_function_uniform;
hfunc2 = @(x) ((1+fN)*((Con-1)*x + 1))^hill/(K^hill+ ((1+fN)*(Con-1)*x + 1)^hill) - x;
for idx=1:3
    fp(idx) = fzero(hfunc2, x0(idx));
end

% homogenization condition
cells_temp = fp(1)*ones(N,1);
cells_temp2 = fp(3)*ones(N,1);
cells_temp(1) = fp(3);
cells_temp2(1) = fp(1);

sigma_c_OFF = std(cells_temp);
sigma_c_ON = std(cells_temp2);

%hin = figure();
%update_cell_figure_continuum_with_I2(hin, pos, a0, cells_temp2, cell_type, 0, 0);
%}
%% Load data and plot

for var=1:numel(a0list)
    a0 = a0list(var);
    Rcell = 0.2*a0;
    
    % file name pattern
    %fpattern = strrep(sprintf('N%d_a0%d_K%d_Con%d_hill%.2f_t%s_xmeanf_%s_%s-v%s',...
    %            N, 10*a0, K, Con, hill, straux, straux2, initialID, straux), '.','p');
    %fpattern = strrep(sprintf('N%d_a0_%.2f_K%d_Con%d_hill%.2f_t%s_xmeanf_%s_%s-v%s',...
    %            N, a0, K, Con, hill, straux, straux2, initialID, straux), '.','p');
    fpattern = strrep(sprintf('N%d_a0_%.2f_K%d_Con%d_hill%.2f_xmeanf_%s_prec_%d_tmax%d%s-v%s',...
                N, a0, K, Con, hill, straux2, prec, tmax, initialID, straux), '.','p');     
    % Initialize all variables and figures
    h1 = figure(1); % <Xi>, sigma(X) trajectories 
    hold on
    h2 = figure(2); % I2 trajectories 
    hold on
    xmean_ini = [];
    xmean_end = [];
    xstd_ini = [];
    xstd_end = [];
    I2_ini = [];
    I2_end = [];
    t_end = [];
    nruns = 0;
    for i = 1:numel(names)
        % first get the filename given by the student
        [tokens, ~] = regexp(names{i},fpattern,'tokens','match');
        if numel(tokens) > 0
            disp(names{i}) % displays the file name
            nruns = nruns + 1;
            % load the data
            load(fullfile(path,strcat(names{i},'.mat')), 'cells_hist', 'fN', 'xmean', 'xstd', 'I', 'Theta', 't');

            I2 = Theta - (2.*xmean - 1).^2;
            xmean_ini(end+1) = xmean(1);
            xmean_end(end+1) = xmean(end);
            xstd_ini(end+1) = xstd(1);
            xstd_end(end+1) = xstd(end);
            I2_ini(end+1) = Theta(1) - (2.*xmean(1) - 1).^2; %I(1)*xstd(1);
            I2_end(end+1) = Theta(end) - (2.*xmean(end) - 1).^2; % I(end)*xstd(end);
            t_end(end+1) = t;

            figure(h1);
            plot(0:t, xmean, 'b-');
            plot(0:t, xstd, 'r-');

            figure(h2);
            plot(0:t, I2, 'b-');        
        end
    end

    % Plot the initial and final points in both figures
    h1 = figure(1);
    figure(h1);

    scatter(t_end, xmean_end, 'bx')
    scatter(t_end, xstd_end, 'rx')
    scatter(zeros(numel(xmean_ini),1), xmean_ini, 'bo')
    scatter(zeros(numel(xstd_ini),1), xstd_ini, 'ro')

    figure(h2)
    scatter(t_end, I2_end, 'bx');
    scatter(zeros(size(I2_ini)), I2_ini, 'bo');

    % Plot the maximum I estimated with first neighbour approx.
    %{
    figure(h1)
    fa0 = sinh(Rcell)*exp(Rcell-a0)/a0;
    Imax = @(p) (6 - 4./sqrt(p*N) - 6./(p*N) - 6*p)*fa0./(1-p)/fN;
    Imax2 = @(p) (6*p - 4./sqrt((1-p)*N) - 6./((1-p)*N))*fa0./p/fN;
    %plot(piv, Imax(piv), 'b')
    %plot(piv, max(Imax(piv), Imax2(piv)), '--b', 'Linewidth', 1.5)
    ylim([-0.05 1])
    %}

    % Set fonts and labels of the map
    figure(h1)
    hold off
    set(gca,'FontSize', 20)
    xlabel('t', 'FontSize', 24)
    ylabel('X', 'FontSize', 24)
    %ylabel(c, 'h=H/N', 'FontSize', 24)
    %ylim([-0.05 0.3])

    % save the map
    qsave = 1;
    if qsave
        name = strrep(sprintf('N%d_a0%.2f_K%d_Con%d_hill%.2f_%s_runs%d_tmax%d',...
                    N, a0, K, Con, hill, initialID, nruns, tmax), '.','p');
        out_file = fullfile(pwd, 'figures', 'time_evolution',... 
            strcat(name,'_x_mean_std')); %filename
        save_figure_pdf(h1, 10, 8, out_file);
        save_figure_eps(h1, 10, 8, out_file);
    end

    % Set fonts and labels of the hamiltonian graph
    figure(h2)
    hold off
    set(gca,'FontSize', 20)
    xlabel('Time (steps)', 'FontSize', 24)
    ylabel('$$I_2$$', 'FontSize', 24)
    %ylim([-0.05 0.3])
    %xlim([0 1500]);
    box on

    % save the pdf
    qsave = 1;
    if qsave
        name = strrep(sprintf('N%d_a0%.2f_K%d_Con%d_hill%.2f_%s_runs%d_tmax%d',...
                    N, a0, K, Con, hill, initialID, nruns, tmax), '.','p');
        out_file = fullfile(pwd, 'figures', 'time_evolution',...
            strcat(name,'_I2')); % filename
        save_figure_pdf(h2, 10, 8, out_file);
        save_figure_eps(h2, 10, 8, out_file);
    end
    %% Plot distribution of final p or I
    %{
    h3 = figure(3);
    idx = (xmean_end~=0); %throw away zero entries
    histogram(xmean_end(idx), 0:0.05:1, 'Normalization', 'count');
    xlim([0 1]);
    xlabel('$$p_{final}$$', 'FontSize', 24);
    ylabel('Count', 'FontSize', 24);
    set(gca,'FontSize', 20);

    %
    h4=figure(4);
    histogram(xstd_end(idx), -0.2:0.05:1, 'Normalization', 'count');
    xlim([-0.2 1]);
    xlabel('$$I_{final}$$', 'FontSize', 24);
    ylabel('Count', 'FontSize', 24);
    set(gca,'FontSize', 20);

    % Save figures
    qsave = 0;
    if qsave
        name = strrep(sprintf('N%d_a0_%.2f_K%d_Con%d_hill%.2f_%s_runs%d',...
                    N, a0, K, Con, hill, initialID, nruns), '.','p');
        out_file = fullfile(pwd, 'figures',...
        'finite_Hill_final_pI_distribution', strcat(name,'_p_distribution'));
        out_file2 = fullfile(pwd, 'figures',...
            'finite_Hill_final_pI_distribution', strcat(name,'_I_distribution'));
        save_figure_pdf(h3, 10, 6, out_file);
        save_figure_pdf(h4, 10, 6, out_file2);
    end
    %}
    %%
    close all;
end