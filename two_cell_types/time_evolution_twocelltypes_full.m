function [Non1_final, Non2_final, I_final, tmax_out, h_decrease, theta_final] = time_evolution_twocelltypes_full(N, dist,...
    a0, Con, K, Rcell, frac_type1, Mcomm, p, tmax, savepath, save_trajectory)
    %{
    % Set parameters of the system
    gridsize = 15;
    N = gridsize^2;
    frac_type1 = 0.5;
    a0 = 0.5;
    Rcell = 0.2*a0;

    % Parameters type 1
    Con1 = 20; 
    K1 = 12;
    p1 = 0.3;
    Rcell1 = 0.5*a0;

    % Parameters type 2
    Con2 = 20; % Con1/2;
    K2 = 8; % K1/2 - K1/5;
    p2 = 0.2;
    Rcell2 = 0.2*a0;
    %}
    %--------------------------------------------------------------------------
    % Check input communication matrix
    %{
    if ~all(size(Mcomm)==[2 2])
        disp('Wrong communication matrix size!');
        return
    end
    for idx=1:numel(Mcomm)
        if ~(Mcomm(idx)==0 || Mcomm(idx)==1)
            disp('Wrong communication matrix entries!');
            return
        end
    end
    %}
    
    % Distance and position
    %[dist, ~] = init_dist_hex(gridsize, gridsize);

    % cell types
    cell_type = zeros(N,1);
    N1 = round(frac_type1*N);
    N2 = N - N1;
    idx1 = randperm(N,N1);
    idx2 = setdiff(1:N, idx1);
    cell_type(idx2) = 1;
    % extra randomization of cell positions by requiring that cell_type also
    % has a Moran's I value of around 0
    [cell_type, ~] = generate_I_new(cell_type, 0, 0.01, dist, a0);
    idx2 = find(cell_type);
    idx1 = setdiff(1:N, idx2);
    
    % Matrix of cell reading
    M = zeros(size(dist));
        % same type
    M(idx1, idx1) = Mcomm(1,1)*sinh(Rcell(1))./(a0*dist(idx1, idx1)).*exp(Rcell(1)-a0*dist(idx1, idx1));
    M(sub2ind(size(M), idx1, idx1)) = Mcomm(1,1)*1; % 1: normal, 0: no self communication

    M(idx2, idx2) = Mcomm(2,2)*sinh(Rcell(2))./(a0*dist(idx2, idx2)).*exp(Rcell(2)-a0*dist(idx2, idx2));
    M(sub2ind(size(M), idx2, idx2)) = Mcomm(2,2)*1; % 1: normal, 0: no self communication
        % different type
    M(idx1, idx2) = Mcomm(1,2)*Rcell(2)/(Rcell(1))*sinh(Rcell(1))./(a0*dist(idx1, idx2)).*exp(Rcell(2)-a0*dist(idx1, idx2)); % conc. cell of type 1 senses due to cells of type 2
    M(idx2, idx1) = Mcomm(2,1)*Rcell(1)/(Rcell(2))*sinh(Rcell(2))./(a0*dist(idx2, idx1)).*exp(Rcell(1)-a0*dist(idx2, idx1));
        % self-communication
    %M(sub2ind(size(M), 1:N, 1:N)) = 1; % 1: normal, 0: no self communication
    %%
    % Interaction strengths
    %f11 = sum(sum(M(idx1, idx1)))/N;
    %f12 = sum(sum(M(idx1, idx2)))/N; 
    %f21 = sum(sum(M(idx2, idx1)))/N; %not equal to f12 if cell radii different
    %f22 = sum(sum(M(idx2, idx2)))/N;

    %dist_vec = a0*dist(1,:);
    %r = dist_vec(dist_vec>0); % exclude self influence
    %fN = sum(sinh(Rcell)*sum(exp(Rcell-r)./r)); % calculate signaling strength
    %%
    % default: tmax not reached
    tmax_out = 0;
    
    % initialize ON cells
    cells = zeros(N,1);
    cells(idx1(randperm(N1,round(p(1)*N1)))) = 1;
    cells(idx2(randperm(N2,round(p(2)*N2)))) = 1;
    
    % initialize figure
    %hin = figure(1);
    t = 0;
    I = [];
    theta11 = [];
    theta12 = [];
    theta22 = [];
    Non1 = [];
    Non2 = [];
    h = [];
    cells_hist = {};
    %update_cell_figure(hin, pos, a0, cells, cell_type, t);
    cells_hist{end+1} = cells;
    Non1(end+1) = sum(cells(idx1));
    Non2(end+1) = sum(cells(idx2));
    I(end+1) = moranI(cells, a0*dist);
    [theta11(end+1), theta12(end+1), theta22(end+1)] = moranI_twocelltypes(cells, idx1, idx2, a0*dist);
    [cells_out, changed, h(end+1)] = update_cells_twocells(cells, idx1, idx2, M, Con, K);
    while changed
        t = t+1;
        if t==tmax
            tmax_out = 1; % tmax reached
            break;
        end
        %k = waitforbuttonpress;
        %update_cell_figure(hin, pos, a0, cells_out, cell_type, t);
        cells_hist{end+1} = cells_out;
        I(end+1) = moranI(cells_out, a0*dist);
        [theta11(end+1), theta12(end+1), theta22(end+1)] = moranI_twocelltypes(cells_out, idx1, idx2, a0*dist);
        Non1(end+1) = sum(cells_out(idx1));
        Non2(end+1) = sum(cells_out(idx2));
        cells = cells_out;
        [cells_out, changed, h(end+1)] = update_cells_twocells(cells, idx1, idx2, M, Con, K);
    end
    
    %% Values to return
    Non1_final = Non1(end);
    Non2_final = Non2(end);
    I_final = I(end); %return value
    theta11_final = theta11(end);
    theta12_final = theta12(end);
    theta22_final = theta22(end);
    theta_final = [theta11_final theta12_final theta22_final];
    h_decrease = all((h(2:end)-h(1:end-1))>0); %whether h decreases uniformly

    %% Save trajectory
    if save_trajectory
        %dir = fullfile(pwd, 'data', 'time_evolution');
        %dir = 'K:\bn\hy\Shared\Yiteng\Multicellularity\data_twocelltypes\time_evolution';
        dir = savepath; %'L:\HY\Shared\Yiteng\Multicellularity\data\time_evolution_twocelltypes';
        fname_str = strrep(...
            sprintf('N1_%d_N2_%d_p1_%.2f_p2_%.2f_a0_%.2f_K1_%d_K2_%d_Con1_%d_Con2_%d_Rcell1_%.2f_Rcell2_%.2f_t%d', ...
            N1, N2, p(1), p(2), a0, K(1), K(2), Con(1), Con(2), Rcell(1), Rcell(2), t), '.', 'p');
        i = 1;
        fname = fullfile(dir, ...
            strcat(fname_str,'-v',int2str(i),'.mat'));
        while exist(fname, 'file') == 2
            i=i+1;
          fname = fullfile(dir, ...
              strcat(fname_str,'-v',int2str(i),'.mat'));
        end
        save(fname, 'N', 'N1', 'a0', 'Con', 'K', 'Rcell', 't', 'Non1',...
            'Non2', 'I', 'h', 'theta11', 'theta12', 'theta22')
    end