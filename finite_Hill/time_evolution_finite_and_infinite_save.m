% Time evolution of a system without noise
% Run and save (1) finite Hill coefficient and (2) infinite Hill
% coefficient with same starting conditions
% NB Probably infinite n not needed, since it's in autonomy for almost all
% interesting cases
close all
clear all
warning off
%%
% Lattice parameters
gridsize = 11;
N = gridsize^2;
%a0 = 6;
%Rcell = 0.2*a0;

% circuit parameters
K = 8;
Con = 16;
hill = 2; % Hill coefficient
prec = 5; %precision

% Initial configuration
initialID = 'binaryrand'; %'binary'; %'uniform' 'allON' 'singleOFF' 'fixedp'
p = 0.5;
sigma = 0.1; 
iniON = round(p*N);

% simulation parameters
tmax = 1000;

% generate cell_type (0 case type 1, 1 case type 2)
cell_type = zeros(N,1); % all the same here
%% Plot phase diagram 
%{
% For single cell in uniform lattice
% Bistability example: K=8, Con=16, a0=1.5, hill=2
set(0,'defaulttextinterpreter', 'latex');
h=figure();
hold on
S = linspace(1+fN,Con*(1+fN),200);
fS = 1 + (Con-1)*S.^hill./(K^hill+S.^hill);
plot(linspace(0, 1, 200), (S-1)/(Con-1), 'b-', 'LineWidth', 1.5);
plot(linspace(0, 1, 200), (fS-1)/(Con-1), 'r-', 'LineWidth', 1.5);
%xlabel('Sensed concentration $$Y_i$$', 'FontSize', 20);
%ylabel('Secreted concentration $$f(Y_i)$$', 'FontSize', 20);
xlabel('$$X_i(t)$$');
ylabel('$$X_i(t+1)$$');
legend({'X_i', 'f(X_i)'}, 'FontSize', 20, 'Location', 'nw');
set(gca, 'FontSize', 24);
xlim([0 1]);
ylim([0 1]);

%Save fig
fname = strrep(sprintf('single_cell_dynamics_Con%d_K%d_hill%.2f_a0_%.2f', Con, K, hill, a0), '.','p');
out_file = fullfile(pwd, 'latex', '17-07-14_finite_Hill', 'figures', fname);
%out_file = fullfile(pwd, 'figures', 'finite_Hill_autonomy_transition', fname);
%save_figure_pdf(h, 10, 8, out_file)
%}
%% repeat many simulations
dx = 0.005;
a0list = 5.05:0.05:5.90; %5.25+dx:dx:5.30-dx;
nruns = 400;
allstable = zeros(numel(a0list), nruns);
for idx = 1:numel(a0list)
    a0 = a0list(idx);
    Rcell = 0.2*a0;
        
    % use hexagonal lattice
    [dist, pos] = init_dist_hex(gridsize, gridsize);
    dist = round(dist, 5);
    dist_vec = a0*dist(1,:);
    r = dist_vec(dist_vec>0); % exclude self influence
    fN = sum(sinh(Rcell)*sum(exp(Rcell-r)./r)); % calculate signaling strength
    gN = sum(sum((sinh(Rcell)*exp(Rcell-r)./r).^2)); % calculate signaling strength

    irun = 1;
    while irun <= nruns
        fprintf('a0 = %.3f, run = %d', a0, irun);
        %% initialize cells
        cells_hist = {};
        if strcmp(initialID, 'uniform')
            cells_hist{1} = rand(N, 1); %uniformly distributed
        elseif strcmp(initialID, 'singleOFF')
            cells_hist{1} = ones(N, 1); 
            x = randi(N);
            cells_hist{1}([x x+1 x+gridsize x+gridsize+1]) = 0;
        elseif strcmp(initialID, 'fixedp')
            cells_hist{1} = p + sigma*randn(N, 1);
            cells_hist{1}(cells_hist{1} < 0) = 0;
            cells_hist{1}(cells_hist{1} > 1) = 1;
        elseif strcmp(initialID, 'fixedp_unif')
            cells_hist{1} = p*ones(N, 1); 
        elseif strcmp(initialID, 'binary')
            cells_temp = zeros(N, 1);
            cells_temp(randperm(N, iniON)) = 1;
            cells_hist{1} = cells_temp;
        elseif strcmp(initialID, 'binaryrand')
            cells_hist{1} = randi(2, N, 1)-1;
        end
        %% Finite Hill simulation
        cells = cells_hist{1};

        % initialize figure
        t = 0;
        %H = [];
        %xmean = [];
        %xstd = [];
        %I = [];
        %xmean(end+1) = mean(cells);
        %xstd(end+1) = std(cells);
        %[I(end+1), ~] = moranI(cells, a0*dist);
        [cells_out, changed, ~] = update_cells_continuum(cells, dist, Con, K, a0, Rcell, hill, prec);

        while changed && t<tmax
            t = t+1;
            cells_hist{end+1} = cells_out;
            %xmean(end+1) = mean(cells_out);
            %xstd(end+1) = std(cells_out);
            %[I(end+1), ~] = moranI(cells, a0*dist);
            cells = cells_out;
            [cells_out, changed, ~] = update_cells_continuum(cells, dist, Con, K, a0, Rcell, hill, prec);
        end
        % Check stability of final state
        [stable, ~] = stability_fun(cells_out, a0, Rcell, K, Con, hill, dist);
        allstable(idx, irun) = stable;
        %% Infinite Hill coefficient simulation
        %{
        %Con2 = fp(3)/fp(1);
        %K2 = K/(fp(1)*Con);
        Con_inf = Con;
        K_inf = K;
        cells = cells_hist{1};
        cells_hist_inf = {}; 
        cells_hist_inf{end+1} = cells;

        % initialize figure
        t = 0;
        %H_inf = [];
        %I_inf = [];
        %[I(end+1), ~] = moranI(cells, a0*dist);
        [cells_out, changed, ~] = update_cells(cells, dist, Con_inf, K_inf, a0, Rcell);
        while changed && t<tmax
            t = t+1;
            cells_hist_inf{end+1} = cells_out;
            %[I_inf(end+1), ~] = moranI(cells, a0*dist);
            cells = cells_out;
            [cells_out, changed, ~] = update_cells(cells, dist, Con_inf, K_inf, a0, Rcell);
        end
        %}
        %% Save whole trajectory if interesting
        qsave = 1;
        if qsave && stable
            irun = irun + 1; % IMPORTANT updating index
            fname_str = strrep(sprintf('N%d_a0_%.2f_Con%.2f_K%.2f_hill%.2f_%s',...
                N, a0, Con, K, hill, initialID), '.','p');
            ifile = 1;
            fname = fullfile(pwd, 'data', 'dynamics_transition 2017-10-25',... 
                strcat(fname_str,'-v',int2str(ifile),'.mat'));
            while exist(fname, 'file') == 2
                ifile = ifile+1;
                fname = fullfile(pwd, 'data', 'dynamics_transition 2017-10-25',...
                    strcat(fname_str,'-v',int2str(ifile),'.mat'));
            end
            cells_ini = cells_hist{1};
            cells_final = cells_hist{end};
            save(fname, 'cells_ini', 'cells_final', 'stable');
        end
    end
end