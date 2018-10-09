% Tests the validity of the EOM for a given trajectory by computing the
% likelihood of the trajectory according to the EOM (transition matrix)
clear variables
close all
warning off
set(0, 'defaulttextinterpreter', 'latex');

%% Load simulation from file
% Load simulated data
%{
folder = 'H:\My Documents\Multicellular automaton\figures\two_signals\macroscopic\period3';
fname_str = 'N225_a0_1p80_rcell_0p20_lambda12_1p20_M_int_0_1_-1_-1_K_0_12_13_8_Con_18_16_iniON_57_56_56_56_EOM';
load(fullfile(folder,fname_str), 'gz', 'N', 'a0', 'rcell', 'Rcell', 'M_int',...
    'Con', 'Coff', 'K', 'lambda', 'hill', 'noise', 'iniON', 'pos', 'dist', 'fN', 'gN');
%}

%folder = 'H:\My Documents\Multicellular automaton\figures\two_signals\macroscopic\chaotic';
folder = 'H:\My Documents\Multicellular automaton\figures\two_signals\cycle_prob\chaotic';
fname_str = 'N225_a0_1p50_rcell_0p20_lambda12_1p00_M_int_0_1_-1_1_K_0_8_11_4_Con_18_16_p_ini_0p5_0p5-v1';
load(fullfile(folder,fname_str)); %, 'cells_hist', 'save_consts_struct', 'positions', 'distances', 't_out');

% parameters
s = save_consts_struct;
N = s.N;
gz = sqrt(N);
a0 = s.a0;
rcell = s.rcell;
Rcell = rcell*a0;
M_int = s.M_int;
Con = s.Con;
Coff = s.Coff; 
K = s.K;
lambda12 = s.lambda12;
lambda = [1 lambda12];
hill = s.hill;
noise = s.noise;

pos = positions;
dist = distances;

% Calculate interaction strength
dist_vec = a0*dist(1,:);
r = dist_vec(dist_vec>0); % exclude self influence
fN = zeros(1, 2);
gN = zeros(1, 2);
for i=1:2
    fN(i) = sum(sinh(Rcell) * sum(exp((Rcell-r)./lambda(i)).*(lambda(i)./r)) ); % calculate signaling strength
    gN(i) = sum(sinh(Rcell)^2 * sum(exp(2*(Rcell-r)./lambda(i)).*(lambda(i)./r).^2 ) ); % calculate noise variance strength
end


%% Calculate transition probabilities

for t=0:t_out-1
    cells_in = cells_hist{t+1};
    cells_out = cells_hist{t+2};
    
    % convert cells into macrostate
    cells_idx = cells_in*[1; 2];
    n_in = histcounts(categorical(cells_idx), {'0', '1', '2', '3'});
    n_in = reshape(n_in, 2, 2);
    
    cells_idx = cells_out*[1; 2];
    n_out = histcounts(categorical(cells_idx), {'0', '1', '2', '3'});
    n_out = reshape(n_out, 2, 2);
    
    % calculate transition probability
    
end

%% Plot mean hamming distance and periodicity probability
h = figure;

plot(1:max_period, mean(d_all, 2)/N, 'go-', 'LineWidth', 1.5, 'MarkerSize', 8);

xlabel('Time');
ylabel('Transition probability');

ylim([0 1]);
set(gca, 'FontSize', 24);
set(h, 'Units', 'Inches', 'Position', [1 1 10 8]);

%% Save
%
% Save data
M_int_s = sprintf('%d_%d_%d_%d', M_int(1,1), M_int(1,2), M_int(2,1), M_int(2,2));
%a0_s = sprintf('%.2f', a0);
%R_s = sprintf('%.2f', rcell);
K_s = sprintf('%d_%d_%d_%d', K(1,1), K(1,2), K(2,1), K(2,2));
Con_s = sprintf('%d_%d', Con(1), Con(2));
%pij_s = sprintf('%.1f_%.1f_%.1f_%.1f', pij(1,1), pij(1,2), pij(2,1), pij(2,2));
iniON_s = sprintf('%d_%d_%d_%d', iniON(1,1), iniON(1,2), iniON(2,1), iniON(2,2));
%I_s = sprintf('%d', mult_dig*round(I, n_dig));

fname_str = strrep(sprintf(...
    'Transition_probs_N%d_a0_%.2f_rcell_%.2f_lambda12_%.2f_M_int_%s_K_%s_Con_%s_iniON_%s_ntrials%d',...
    N, a0, rcell, lambda(2), M_int_s, K_s, Con_s, iniON_s, ntrials), '.', 'p');
folder = 'K:\bn\hy\Shared\Yiteng\Multicellularity\videos\selected\transition_probs_EOM';

qsave = 0;
if qsave
    save(fullfile(folder, fname_str));
end

% Save figure
qsave = 0;
save_figure(h, 10, 8, fullfile(folder, fname_str), '.pdf', qsave);
%}
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