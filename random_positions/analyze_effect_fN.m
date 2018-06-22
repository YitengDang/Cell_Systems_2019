%% Time evolution of lattice with visualization
clear all
close all
set(0,'defaulttextinterpreter', 'latex');
%%
% geometric parameters
Lx = 1;
n = 12; % nmax = L/R (square packing)
N = n^2; % total number of particles
%N = round(eta*L^2/(pi*R^2));
Ly = sqrt(3)/2*Lx;
r_av = sqrt(Lx*Ly/N)/2; % estimate for NND of random distribution (Clark & Evans, 1954)
a0 = Lx/n;
rcell = 0.2;
R = rcell*a0; % disc radius
eta = N*pi*R^2/(Lx*Ly); %packing fraction 
if eta > pi/(2*sqrt(3)) % max. packing fraction
    disp('packing fraction too high! Abort evaluation');
    pause(10);
end

% random lattice (reference)
[pos_rand, dist_rand] = initial_cells_random_periodic_alt(n, Lx, Ly, R); % random config
    
% loop over mcsteps
mcsteps_all = 10.^(0:6);
fN_all = zeros( numel(mcsteps_all), N );
for idx_mc = 1:numel(mcsteps_all)
    mcsteps = mcsteps_all(idx_mc);
    
    % Initial configuration
    cells = ones(N, 1);

    [pos, dist, fN0] = initial_cells_random_markov_periodic(n, Lx, R, mcsteps);

    %% draw figure
    %{
    t=0;
    hin = figure(1);
    update_figure_periodic_long(pos+R, N, Lx, Ly, R, cells, t)

    h2 = figure(2);
    update_figure_periodic_long(pos_rand+R, N, Lx, Ly, R, cells, t)
    %}
    %% Check that distances are the same
    % plot distance distribution
    %{
    uniq = unique(round(dist(1,:), 10));
    figure;
    hold on
    for i=1:N
        %uniq = unique(dist(i,:));
        %idx = dist(i,:)==uniq(2);
        %plot(0, dist(i, idx));
        %plot(i, dist(i,:), 'bo');
        disp(all(uniq == unique(round(dist(i,:), 10)) ) )
    end
    %}

    %% Nearest neighbour distances
    %{
    dist2 = reshape(dist(dist>0), [N-1 N]);
    nnd = min(dist2);
    nnd_av = mean(nnd);

    % Histogram of NND
    h=figure();
    hold on
    nbins = 20;
    histogram(nnd, linspace(2*R,(max(nnd)-2*R), nbins), 'Normalization', 'pdf');
    xlabel('NND');
    ylabel('Probability');
    title(sprintf('$$ \\langle NND \\rangle = %.2f$$, target $$= %.2f $$', nnd_av, r_av));
    set(gca, 'FontSize', 24);
    set(h, 'Position', [500 200 840 720]);

    % target a0
    plot([r_av r_av], [0 100], 'r--');
    % minimum distance
    plot([2*R 2*R], [0 100], 'k--');

    % Clark & Evans calculation
    %{
    rho = @(r) 2*pi*N/(Lx*Ly)*r.*exp(-pi*N/(Lx*Ly)*r.^2);
    rvals = linspace(0, max(nnd)*1.2, 1000);
    plot(rvals, rho(rvals), 'LineWidth', 2);
    %xlim([0 1]);
    %}
    %}
    %% Plot fN distribution
    fN = zeros(N, 1);
    fN_rand = zeros(N, 1);

    for i=1:N
        dist_vec = dist(i,:)*(Lx/n);
        r = dist_vec(dist_vec>0); % exclude self influence
        fN(i) = sum(sinh(R)*sum(exp(R-r)./r)); % calculate signaling strength

        dist_vec = dist_rand(i,:)*(Lx/n);
        r = dist_vec(dist_vec>0); % exclude self influence
        fN_rand(i) = sum(sinh(R)*sum(exp(R-r)./r)); % calculate signaling strength
    end
    %{
    figure;
    hold on
    plot(1:N, fN/fN0, 'o-');
    plot(1:N, fN_rand/fN0, 'o-');
    plot([1 N], [mean(fN)/fN0 mean(fN)/fN0], '--');
    plot([1 N], [mean(fN_rand)/fN0 mean(fN_rand)/fN0], '--');
    title({sprintf('MC: $$\\sigma(f_N) = %.5f$$', std(fN)),...
        sprintf('Brute force: $$\\sigma(f_N) = %.5f$$', std(fN_rand))})
    ylim([0 1.1]);
    ylabel('fN/fN0');
    legend({'MC', 'Brute force'});
    %}
    h=figure(3);
    hold on
    %nbins = 10;
    edges = linspace(0.8, 1.25, 15);
    map = colormap();
    histogram(fN/fN0, edges, 'FaceColor', map(1,:), 'Normalization', 'pdf');
    histogram(fN_rand/fN0, edges, 'FaceColor', map(end,:),'Normalization', 'pdf');
    plot([mean(fN)/fN0 mean(fN)/fN0], [0 40], '--', 'Color', map(1,:));
    plot([mean(fN_rand)/fN0 mean(fN_rand)/fN0], [0 40], 'r--', 'Color', map(end,:));
    legend({sprintf('MC 10^{%d} steps', log10(mcsteps)), 'Brute force'});
    xlabel('$$f_N/f_N^0$$');
    ylabel('Probability');
    set(gca,'FontSize', 24);
    %ylim([0 7]);
    %% Save figure
    fname_str = strrep(sprintf('fN_distribution_compare_N%d_rcell%.2f_hist_mcsteps%d',...
        N, rcell, mcsteps), '.', 'p');
    folder = 'H:\My Documents\Multicellular automaton\figures\random_positions';
    fname = fullfile(folder, fname_str);
    qsave = 1;
    if qsave
        save_figure(h, 10, 8, fname, '.pdf');
    end
    %%
    fN_all(idx_mc, :) = fN;
    close all
end
%% Save data

folder = 'H:\My Documents\Multicellular automaton\data\random_positions\fN_effect';
fname_str = strrep(sprintf('fN_vs_mcsteps_N%d_rcell%.2f',...
    N, rcell, mcsteps), '.', 'p');
fname = fullfile(folder, fname_str);
save(fname, 'mcsteps_all', 'fN_all', 'Lx', 'n', 'rcell');

%% Plot fN vs mcsteps
h4 = figure(4);
errorbar(log10(mcsteps_all), mean(fN_all/fN0, 2), std(fN_all/fN0, 1, 2), 'LineWidth', 2);
xlabel('$$log_{10}$$(MC steps)');
ylabel('$$f_N/f_N^0$$');
set(gca,'FontSize', 24);

folder = 'H:\My Documents\Multicellular automaton\figures\random_positions';
fname_str = strrep(sprintf('fN_vs_mcsteps_N%d_rcell%.2f',...
    N, rcell, mcsteps), '.', 'p');
fname = fullfile(folder, fname_str);
qsave = 1;
if qsave
    save_figure(h4, 10, 8, fname, '.pdf');
end