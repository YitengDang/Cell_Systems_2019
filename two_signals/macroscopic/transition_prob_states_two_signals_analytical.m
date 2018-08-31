% Calculates the state transition probabilities for two signals through
% analytical theory
% pij: for transitions between (i,j) states
clear variables
close all
warning off
set(0, 'defaulttextinterpreter', 'latex');
%%
% Parameters of the system
gz = 15;
N = gz^2;
a0 = 1.5;
rcell = 0.2;
Rcell = rcell*a0;

M_int = [0 1; -1 -1];
Con = [18 16];
Coff = [1 1];
K = [0 12; 13 8];
lambda = [1 1.2];
hill = Inf;
noise = 0;
%}

% Load from file
%{
%folder = 'H:\My Documents\Multicellular automaton\app\data\time_evolution\sample_trajectories\typeII';
folder = 'H:\My Documents\Multicellular automaton\data\two_signals\transition_matrix\sim';
fname_str = 'W_sim_M_int_0_1_-1_-1_N225_a0_1p50_rcell_0p20_K_0_12_13_8_Con_18_16_iniON_50_61_58_56_delta_t_1_trials_1000';
load(fullfile(folder,fname_str), 'gz', 'N', 'a0', 'rcell', 'Rcell', 'M_int',...
    'Con', 'Coff', 'K', 'lambda', 'hill', 'noise', 'iniON');
%}

% (1) original hexagonal lattice
%[pos,ex,ey] = init_cellpos_hex(gridsize,gridsize);
%dist = dist_mat(pos,gridsize,gridsize,ex,ey);
% (2) new lattice
mcsteps = 0;
[pos, dist] = initial_cells_random_markov_periodic(gz, mcsteps, rcell);

% Calculate interaction strength
dist_vec = a0*dist(1,:);
r = dist_vec(dist_vec>0); % exclude self influence
fN = zeros(1, 2);
gN = zeros(1, 2);
for i=1:2
    fN(i) = sum(sinh(Rcell) * sum(exp((Rcell-r)./lambda(i)).*(lambda(i)./r)) ); % calculate signaling strength
    gN(i) = sum(sinh(Rcell)^2 * sum(exp(2*(Rcell-r)./lambda(i)).*(lambda(i)./r).^2 ) ); % calculate noise variance strength
end

% background p^{ij}
%pij = [0.2 0.2; 0.3 0.3]; % (0,0), (0,1); (1,0), (1,1)
iniON = [220 5; 0 0];
pij = iniON/N;
if sum(pij(:))~=1
    warning('probabilities do not add up to 1!');
end
%% Calculate W
W = transition_prob_states_two_signals_analytical_calc(M_int, Con, Coff, K, fN, gN, pij);

%% Plot matrix
h = figure(1);
imagesc(1:4, 1:4, W')
hold on
for i=1:3
    line([0.5, 4.5], [i+0.5, i+0.5], 'Color', 'r', 'LineWidth', 1);
    line([i+0.5, i+0.5], [0.5, 4.5], 'Color', 'r', 'LineWidth', 1);
end
c = colorbar;
set(gca, 'Ydir', 'normal', 'FontSize', 20)
xlabel('$$X(t)$$', 'FontSize', 24)
ylabel('$$X(t+1)$$', 'FontSize', 24)
c.Label.String = 'Probability';
set(gca, 'XTick', 1:4, 'YTick', 1:4, 'XTickLabels', ...
    {'(0,0)','(1,0)','(0,1)','(1,1)'},...
    'YTickLabels', {'(0,0)','(1,0)','(0,1)','(1,1)'});
set(h, 'Position', [100 100 640 520]);
caxis([0 1]);
%}

%% Save
% Save data
folder = 'H:\My Documents\Multicellular automaton\data\two_signals\transition_matrix';
M_int_s = sprintf('%d_%d_%d_%d', M_int(1,1), M_int(1,2), M_int(2,1), M_int(2,2));
a0_s = sprintf('%.2f', a0);
R_s = sprintf('%.2f', rcell);
K_s = sprintf('%d_%d_%d_%d', K(1,1), K(1,2), K(2,1), K(2,2));
Con_s = sprintf('%d_%d', Con(1), Con(2));
%pij_s = sprintf('%.1f_%.1f_%.1f_%.1f', pij(1,1), pij(1,2), pij(2,1), pij(2,2));
iniON_s = sprintf('%d_%d_%d_%d', iniON(1,1), iniON(1,2), iniON(2,1), iniON(2,2));
%I_s = sprintf('%d', mult_dig*round(I, n_dig));

%fname_str = strrep(sprintf('W_analytical_M_int_%s_N%d_a0_%s_rcell_%s_K_%s_Con_%s_pij_%s',...
%    M_int_s, N, a0_s, R_s, K_s, Con_s, pij_s), '.', 'p');
fname_str = strrep(sprintf('W_analytical_M_int_%s_N%d_a0_%s_rcell_%s_K_%s_Con_%s_iniON_%s',...
    M_int_s, N, a0_s, R_s, K_s, Con_s, iniON_s), '.', 'p');

qsave = 0;
if qsave
    save(fullfile(folder, 'analytical', fname_str));
end

% Save figure
qsave = 0;
save_figure(h, 10, 8, fullfile(folder, 'figures', fname_str), '.pdf', qsave);

%% 
%{
%function W = transition_prob_states_two_signals_analytical_func(M_int, Con, Coff, K, fN, gN, pij)
    % Calculate transition matrix
    % pij: Background p^{ij}
    W = zeros(4, 4); %transition matrix
    I_in = [0 0];
    % loop over p_in, p_out
    for i=2 %1:4
        [i1, i2] = ind2sub([2 2], i);
        p_in = [i2-1 i1-1]; % order: (0,0), (0,1), (1,0), (1,1)
        %disp(p_in);
        for j=2 %1:4
            [i1, i2] = ind2sub([2 2], j);
            p_out = [i2-1 i1-1];
            
            p_nei = [sum(pij(2, :)) sum(pij(:, 2))];
            %p_in = [0 0];
            %p_out = [0 0];
            % ------ calculate W1 (prob. p_out^{s} = 1)------
            %W1 = zeros(1, 2);

            % self-contributions to Y
            Y_self = Con.*p_in + Coff.*(1-p_in);

            % neighbour contributions to Y
            Y_nei_avg = fN.*(p_nei.*Con + (1-p_nei).*Coff + (Con-1).*(1-p_nei).*I_in); % neighbour contributions

            % total sensed concentration
            Y = Y_self + Y_nei_avg;
            Y_mat = repmat(Y, 2, 1);

            % std
            sigma = sqrt(p_nei.*(1-p_nei).*gN).*(Con-1);
            sigma_mat = repmat(sigma, 2, 1);

            % Note that W1 = Ponon
            %disp( M_int.*(Y_mat - K) );
            %disp(sigma_mat);
            W1 = prod(1 - abs(M_int).*normcdf(0, M_int.*(Y_mat - K), sigma_mat), 2 )';
            % ------ ------------ ------

            P1 = (1-W1).*(1-p_out) + W1.*p_out;
            W(i,j) = prod(P1);
        end
    end

    % Check whether the probabilities add up to 1
    if ~all(sum(W, 2)==1)
        warning('Probabilities do not add up correctly!');
        disp('Probabilities add up to:');
        disp(all(sum(W, 2)==1));
    end
%end
%}