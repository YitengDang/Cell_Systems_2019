% Predicts whether a travelling wave can propagate according to nearest
% neighbour interactions
%{
clear all
close all
set(0,'defaulttextinterpreter', 'latex')

%% Required input
% Manual input
gz = 15;
N = gz^2;
a0 = 1.5;
rcell = 0.2;
Rcell = rcell*a0;
lambda = [1 1.2];
M_int = [1 1; -1 0]; % network 19 
%{
M_int = [0 1; -1 1]; % network 15 reversed
M_int = [1 -1; 1 0]; % network 15
M_int = [1 1; -1 0]; % network 19
M_int = [1 -1; 1 1]; % network 33
M_int = [-1 -1; 1 1]; % network 34
M_int = [-1 1; -1 1]; % network 36
%}
Con = [18 16];
%K = zeros(2);
K = [0 9; 11 4];

% get pos, dist
mcsteps = 0;
[pos, dist] = initial_cells_random_markov_periodic(gz, mcsteps, rcell);
%%
trav_wave_conds_met = temp(gz, a0, dist, rcell, lambda, M_int, Con, K,...
    wave_type, states_perm, num_waves, bandwidth);
%}
function trav_wave_conds_met = travelling_wave_stability_predictor_general_func(gz, a0, dist, rcell, lambda, M_int, Con, K,...
    wave_type, states_perm, num_waves, bandwidth)
    % calculate fN
    Rcell = a0*rcell;
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

    %% Calculate Y_all
    default_states = [0 0; 0 1; 1 0; 1 1];

    % calculate Y(alpha)
    % neighbour data
    switch wave_type
        case 1
            n_nei = [2	0	0	4; % EF
                2	2	0	2; % F
                2	2	2	0; % M 
                0	2	2	2; % B
                0	0	2	4; % EB
                0	0	0	6]; % E
        case 2
            n_nei = [3	0	0	3;
                2	3	0	1;
                1	2	3	0;
                0	1	2	3;
                0	0	1	5;
                0	0	0	6];
        case 3
            n_nei = [1	0	0	5;
                2	1	0	3;
                3	2	1	0;
                0	3	2	1;
                0	0	3	3;
                0	0	0	6];
    end
    states = default_states(states_perm, :); % F, M, B, E
    types = [states(4,:); states(1,:); states(2,:); states(3,:); states(4,:); states(4,:)];
    Y_self = (Con-1).*types + 1;

    % calculate Y_nei
    Y_nei = zeros(size(n_nei, 1), 2);
    for i=1:6
        Y_nei(i,1) = fnn(1)*n_nei(i,:)*((Con(1)-1)*states(:,1)+1);
        Y_nei(i,2) = fnn(2)*n_nei(i,:)*((Con(2)-1)*states(:,2)+1);
    end
    %
    % estimate p
    tmp = num_waves*bandwidth;
    tmp2 = [tmp tmp tmp gz-3*tmp];
    p = (tmp2*states)/gz;

    % calculate Y_mf (mean-field contribution)
    z = 6; % coordination number
    Y_mf = zeros(1, 2);
    Y_mf(1) = (fN(1) - z*fnn(1))*( Con(1)*p(1) + (1-p(1)) );
    Y_mf(2) = (fN(2) - z*fnn(2))*( Con(2)*p(2) + (1-p(2)) );
    %disp(Y_mf)

    Y_all = Y_self + Y_nei + Y_mf;

    %% Check conditions for specific parameter set
    conditions_met = zeros(6,1);
    % (1) E_F->F
    % (2) F->M
    % (3) M->B
    % (4) B->E
    % (5) E_B->E
    % (6) E->E 
    targets = [1 2 3 4 4 4]; % recall F, M, B, E
    output_state_list = zeros(6, 2);
    for i=1:6
        %i = 1;
        Y = Y_all(i, :); % 2x1
        output_state = prod(((Y - K).*M_int>0) + (1-abs(M_int)), 2)';
        target_state = states(targets(i),:);
        conditions_met(i) = all(output_state==target_state);
        output_state_list(i,:) = output_state;
    end
    %disp(output_state_list);
    %disp(conditions_met);

    trav_wave_conds_met = all(conditions_met);
end