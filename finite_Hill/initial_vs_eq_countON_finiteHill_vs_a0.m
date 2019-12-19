% This script calculates the map of p_in p_eq using exact simulation.
close all
clear variables
%warning off
%%
% Parameters of the system
gridsize = 15;
N = gridsize^2;
K = 8;
Con = 16;
hill = 2;
prec = 8; 
noise = 0;
n_smpl = 100;
a0_all = [4 5.2 5.6 10];

% method for sampling initial states; binary=ON/OFF, normal = N(p, sigma)
sampmethod = 'in_out_binary_v2'; %'montecarlo'; %'binary' 'normal' 'uniform'
subfolder = sampmethod; %subfolder for storing graphs
data_subfolder = 'data_binary'; %subfolder for storing data

% use hexagonal lattice
[dist, pos] = init_dist_hex(gridsize, gridsize);

% Calculate the signaling strength
%dist_vec = dist(1, :);
%r = a0*dist_vec(dist_vec>0); % exclude self influence
%fN = sum(Rcell*sum(exp(Rcell-r)./r)); % calculate signaling strength
%kon = N*(K-fN-Con)/(Con-1)/fN;
%koff = N*(K-fN-1)/(Con-1)/fN;

save_fig_folder = 'H:\My Documents\Multicellular automaton\figures\finite_Hill\pin_pout\in_out_binary_v2\vs_a0';
%%
for run=1:numel(a0_all)
    close all;
    
    % loop-dependent variables
    a0 = a0_all(run);
    Rcell = a0*0.2;
    fprintf('a0 = %.2f \n', a0);
    %{
    % calculate signaling strength
    dist_vec = dist(1, :);
    r = a0*dist_vec(dist_vec>0); % exclude self influence
    fN = sum(Rcell*sum(exp(Rcell-r)./r)); 

    % Find uniform lattice fixed points
    fp = zeros(3, 1);
    x0 = [0.03 0.2 0.65]; %estimates based on previous graph
    %hfunc = @update_function_uniform;
    hfunc2 = @(x) ((1+fN)*((Con-1)*x + 1))^hill/(K^hill+ (((1+fN)*(Con-1)*x + 1))^hill) - x;
    for idx=1:3
        fp(idx) = fzero(hfunc2, x0(idx));
    end

    fname_str = strrep(sprintf('pin_pout_N%d_Con%d_K%d_a0_%.2f_hill%.2f_noise%.2f_prec%d_%s_nsmpl%d',...
        gridsize^2, Con, K, a0, hill, noise, prec, sampmethod, n_smpl), '.', 'p');

    % calculate the map
    if strcmp(sampmethod, 'montecarlo')
        [count, t_av, I_av, ntmax] = count_eq_parallel_finiteHill_monte_carlo(...
            dist, Con, K, Rcell, a0, hill, prec);    
    elseif strcmp(sampmethod, 'binary')
        [count, t_av, I_av] = count_eq_parallel_finiteHill_binary(...
            dist, Con, K, Rcell, a0, hill, prec); 
    elseif strcmp(sampmethod, 'in_out_binary')
        [count, t_av, I_av] = count_eq_parallel_finiteHill_fin_fout(...
            dist, Con, K, Rcell, a0, hill, noise, prec, fp, n_smpl); %NB need fixed points
    elseif strcmp(sampmethod, 'in_out_binary_v2')
        [count, t_av, I_av] = count_eq_parallel_finiteHill_fin_fout_v2(...
            dist, Con, K, Rcell, a0, hill, noise, prec, fp, n_smpl); %NB need fixed points
    elseif strcmp(sampmethod, 'normal') 
        [count, t_av, I_av] = count_eq_parallel_finiteHill_normal(...
            dist, Con, K, a0, Rcell, hill);
    elseif strcmp(sampmethod, 'uniform')
        [count, t_av, I_av] = count_eq_parallel_finiteHill_uniform(...
        dist, Con, K, a0, Rcell, hill);
    else 
        disp('ERROR! sampling method not found');
    end
    p = (0:N)./N;
    prob = transpose(count./repmat(sum(count,2),1,N+1));
    %}
    %% save mat file
    %{
    save_folder_data = 'H:\My Documents\Multicellular automaton\figures\finite_Hill\pin_pout\data_binary_v2\vs_a0';
    %file_out = fullfile(pwd, 'figures', 'pin_pout', data_subfolder, strcat(fname_str,'.mat'));
    file_out = fullfile(save_folder_data, fname_str);
    save(file_out);
    %}
    %% Load from file
    %
    fname_str = strrep(sprintf('pin_pout_N%d_Con%d_K%d_a0_%.2f_hill%.2f_noise%.2f_prec%d_%s_nsmpl%d',...
        gridsize^2, Con, K, a0, hill, noise, prec, sampmethod, n_smpl), '.', 'p');
    save_folder_data = 'H:\My Documents\Multicellular automaton\figures\finite_Hill\pin_pout\data_binary_v2\vs_a0';
    load_file = fullfile(save_folder_data, fname_str);
    load(load_file);
    %}
    %% Calculate I
    % Compute mutual information
    % Compute S[ P(p_out) ]
    prob_pout = sum(prob, 2)/(N+1); % P(p_out)
    idx = prob_pout>0;
    S_pout = - sum(prob_pout(idx).*log2(prob_pout(idx)));

    % Compute S[ P(p_out | p_in) ]
    P_pin = 1/(N+1);
    S_pout_pin = zeros(N+1, 1);
    for k=1:N
        P_pout_pin = prob(:, k);
        idx2 = P_pout_pin > 0;
        S_pout_pin(k) = - sum(P_pout_pin(idx2).*log2(P_pout_pin(idx2)));
    end

    % Compute I(pin, pout)
    thisI = S_pout - sum(P_pin.*S_pout_pin);
    
    %% plot the map
    set(0,'defaulttextinterpreter', 'tex');
    h1 = figure(1);
    im_fig = imagesc(p,p,prob);
    % set title and font
    %title(sprintf('$$N = %d, K = %d, S_{ON} = %.1f, a_0 = %.1f, R = %.1f$$', ...
    %    N, K, Con, a0, Rcell),'FontSize', 18)
    title(sprintf('I = %.3f', thisI));
    set(gca,'Ydir','normal','FontSize', 32);
    % set invisible parts where count is zero
    set(im_fig, 'AlphaData', count' > 0);
    % set colorbar and labels
    c = colorbar;
    colormap('winter');
    c.Label.String = 'Probability';
    xlabel('f_{in}')
    ylabel('f_{out}')
    ylim([-0.02 1.02]);
    
    % Organize and save
    save_fig = 1; % save figure? 0: no, 1: yes
    if save_fig > 0
        %out_file = fullfile(pwd, 'figures', 'pin_pout', subfolder, fname);
        out_file = fullfile(save_fig_folder, strcat(fname_str, '_fin_fout_v3') );
        save_figure(h1, 10, 8, out_file, '.pdf', save_fig);
    end

    % Solve for uniform lattice fixed points and include in plot
    %{
    fp = zeros(3, 1);
    x0 = [0.03 0.2 0.65]; %estimates based on previous graph
    hfunc = @update_function_uniform;
    hfunc2 = @(x) hfunc(x, hill, Con, K, fN);
    for idx=1:3
        fp(idx) = fzero(hfunc2, x0(idx));
    end
    figure(h1);
    hold on
    plot([0 1], [fp(1) fp(1)], 'r--'); % stable fixed point
    %plot([0 1], [fp(2) fp(2)], 'r--');
    plot([0 1], [fp(3) fp(3)], 'r--'); % stable fixed point
    plot([fp(2) fp(2)], [0 1], 'g--'); % unstable fixed point

    % Organize and save
    save_fig = 1; % save figure? 0: no, 1: yes
    if save_fig > 0
        out_file = fullfile(pwd, 'figures', 'pin_pout', subfolder, strcat(fname, '_v2') );
        save_figure_pdf(h1, 10, 8, out_file);
        save_figure_eps(h1, 10, 8, out_file);
    end
    %}
    %%
    % Plot the average number of steps it takes to reach equilibrium
    h2 = figure(2);
    plot(p, t_av, 'r-o')
    set(gca,'FontSize', 32)
    %title(sprintf('$$N = %d, K = %d, S_{ON} = %d, a0 = %.1f, R = %.1f$$', ...
    %    N, K, Con, a0, Rcell),'FontSize', 18)
    xlabel('f_{in}', 'FontSize', 32)
    %ylabel('Average # steps for eq.', 'FontSize', 24)
    ylabel('\langle t_{eq} \rangle');
    pmax = p(t_av == max(t_av));
    %title(sprintf('max $$t_{eq}$$ at $$p_{in}$$ = %.3f', pmax));

    save_fig = 0; % save figure? 0: no, 1: yes
    if save_fig > 0
        %out_file = fullfile(pwd, 'figures', 'pin_pout', subfolder, strcat(fname_str, '_tav'));
        out_file = fullfile(save_fig_folder, strcat(fname_str, '_t_av') );
        save_figure(h2, 10, 8, out_file, '.pdf', save_fig);
    end
    %%
    %{
    % Plot the average I in equilibrium
    h3 = figure(3);
    plot(p, I_av, 'r-o')
    set(gca,'FontSize', 32)
    %title(sprintf('$$N = %d, K = %d, S_{ON} = %d, a0 = %.1f, R = %.1f$$', ...
    %    N, K, Con, a0, Rcell),'FontSize', 18)
    xlabel('f_{in}', 'FontSize', 32)
    ylabel('Average I in equilibrium', 'FontSize', 32)
    
    save_fig = 0; % save figure? 0: no, 1: yes
    if save_fig > 0
        %out_file = fullfile(pwd, 'figures', 'pin_pout', subfolder, strcat(fname_str, '_Iav'));
        out_file = fullfile(save_fig_folder, strcat(fname_str, '_I_av') );
        save_figure(h3, 10, 8, out_file, '.pdf', save_fig);
    end
    %}
    %%
    %{
    h4 = figure(4);
    plot(p, ntmax, 'r-o')
    set(gca,'FontSize', 24)
    title(sprintf('$$N = %d, K = %d, S_{ON} = %d, a0 = %.1f, R = %.1f$$', ...
        N, K, Con, a0, Rcell),'FontSize', 18)
    xlabel('$$p_{in}$$', 'FontSize', 24)
    ylabel('Average I in equilibrium', 'FontSize', 24)
    %}
    %{
    h3 = figure(3);
    % calculate the analytical formula and plot
    [~,omegak] = entropy_eq_sphere(dist_vec, Son, K, a0, Rcell);
    plot((0:N)/N , log(omegak), 'LineWidth', 1.5)
    %}

end