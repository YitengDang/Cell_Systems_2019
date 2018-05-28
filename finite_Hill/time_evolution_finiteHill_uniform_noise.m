close all
clear variables

% Parameters of the system
gridsize = 15;
N = gridsize^2;
a0 = 5.6;
Rcell = 0.2*a0;
K = 8;
Con = 16;
hill = 2;
%noise = 10^(0.5);
tmax = 1000;
nsim = 10^4; %number of simulations

% use hexagonal lattice
[dist, pos] = init_dist_hex(gridsize, gridsize);

% Calculate the signaling strength
dist_vec = dist(1, :);
r = a0*dist_vec(dist_vec>0); % exclude self influence
fN = sum(Rcell*sum(exp(Rcell-r)./r)); % calculate signaling strength

% Find zeros
fp = zeros(3, 1);
x0 = [0.03 0.2 0.65]; %estimates based on previous graph
hfunc = @update_function_uniform;
hfunc2 = @(x) hfunc(x, hill, Con, K, fN);
for idx=1:3
    fp(idx) = fzero(hfunc2, x0(idx));
end

% Plot
%{
h1=figure(1);
hold on
x_all = linspace(0,1);
plot(x_all, x_all);
plot(x_all,  ((1+fN)*(Con-1).*x_all + 1).^hill./( ((1+fN)*(Con-1).*x_all + 1).^hill + K^hill )  )
%plot([fp(1) fp(1)], [0 1]);
%plot([fp(2) fp(2)], [0 1]);
%plot([fp(3) fp(3)], [0 1]);
%}

% Minimum noise
%dK = ( (1-fp(2))/fp(2) )^(1/2) * ((1+fN)*(Con-1)*fp(1) + 1) - K;
noise_all = 10.^(-1:0.1:1);
for trial=1:numel(noise_all)
    noise = noise_all(trial);
    x_out = zeros(nsim, tmax+1);
    %x_out(:, 1) = fp(1); %initial state
    x_out(:, 1) = fp(3); %initial state
    dK = noise*randn(nsim, tmax);
    
    %sim = 1;
    %x_state = zeros(nsim, 1);
    for sim=1:nsim
        for t=1:tmax
            this_x_out = (((1+fN)*(Con-1)*x_out(sim, t) + 1))^hill/((K+dK(sim,t))^hill+(((1+fN)*(Con-1)*x_out(sim, t) + 1))^hill);
            x_out(sim, t+1) = this_x_out;
            %switch x_state(sim)
            %    case 0
            %        x_state(sim) = this_x_out > fp(2);
            %    case 1
            %        x_state(sim) = 1 + (this_x_out < fp(2));
            %end
        end
    end
    x_out_final = x_out(:, tmax+1);
    x_out_max = max(x_out, [], 2);
    x_out_min = min(x_out, [], 2);
    x_cross1 = x_out_max > fp(2);
    x_cross2 = x_out_min < fp(2); %zeros(nsim, 1);
    
    %{
    for i=1:nsim
        idx = find(x_out(i, :) > fp(2), 1); 
        if idx
            disp(idx);
            x_cross2(i) = min(x_out(idx+1:end) ) < fp(2);
            disp(min(x_out(idx+1:end)));
        end
    end
    %}
    
    % cross 1
    h3=figure(3);
    hold on
    edges = [linspace(0,fp(1), round(fp(1)/0.01) ) linspace(fp(1), fp(2), round((fp(2)-fp(1))/0.01))...
        linspace(fp(2), fp(3), round((fp(3) - fp(2))/0.01)) linspace(fp(3), 1, round((1-fp(3))/0.01))];
    histogram(x_out_max, edges, 'Normalization', 'pdf');
    plot([fp(1) fp(1)], [0 100]);
    plot([fp(2) fp(2)], [0 100]);
    plot([fp(3) fp(3)], [0 100]);
    title(sprintf('max, frac: %.4f', mean(x_cross1)));
    xlim([0 1]);
    ylim([0 30]);
    set(gca, 'FontSize', 24);
    
    % cross 2
    h4=figure(4);
    hold on
    edges = [linspace(0,fp(1), round(fp(1)/0.01) ) linspace(fp(1), fp(2), round((fp(2)-fp(1))/0.01))...
        linspace(fp(2), fp(3), round((fp(3) - fp(2))/0.01)) linspace(fp(3), 1, round((1-fp(3))/0.01))];
    histogram(x_out_min, edges, 'Normalization', 'pdf');
    plot([fp(1) fp(1)], [0 100]);
    plot([fp(2) fp(2)], [0 100]);
    plot([fp(3) fp(3)], [0 100]);
    title(sprintf('max, frac: %.4f', mean(x_cross2)));
    xlim([0 1]);
    ylim([0 30]);
    set(gca, 'FontSize', 24);
    
    % Print results
    %fprintf('Fraction crossed 1: %.4f \n', mean(x_cross1));
    fprintf('Fraction crossed 2: %.4f \n', mean(x_cross2));
    %% Save data 
    path = 'H:\My Documents\Multicellular automaton\temp';
    fname_str = strrep(sprintf('uniform_lattice_trajectories_N%d_Con%d_K%d_a0_%.2f_hill%.2f_noise_%.3f_tmax%d_nsim%d_cross2', ...
        N, Con, K, a0, hill, noise, tmax, nsim), '.', 'p');
    
    % figures
    save_fig = 0;
    if save_fig > 0
        out_file = fullfile(path, strcat(fname_str, '_final_x_hist'));
        save_figure_pdf(h3, 10, 8, out_file);
        save_figure_eps(h3, 10, 8, out_file);
    end
    
    % data
    close all
    file_out = fullfile(path, fname_str);
    save(file_out, 'a0', 'Con', 'edges', 'fN', 'fp', 'gridsize', 'hill', 'K',...
        'N', 'noise', 'nsim', 'x_cross2', 'x_out_min', 'tmax');
end

%% Visualize single trajectories
%{
h2=figure(2);
hold on
idx = 10;
plot(0:tmax, x_out(idx,:) );
plot([0 tmax], [fp(1) fp(1)]);
plot([0 tmax], [fp(2) fp(2)]);
plot([0 tmax], [fp(3) fp(3)]);
%}
%% Load multiple and plot fraction crossed against noise level
%
path = 'H:\My Documents\Multicellular automaton\temp';
noise_all = 10.^[-1:0.1:1];
frac_final = zeros(numel(noise_all), 1);
frac_cross1 = zeros(numel(noise_all), 1);
frac_cross2 = zeros(numel(noise_all), 1);

for i=1:numel(noise_all)
    noise = noise_all(i);
    fname_in = strrep(sprintf('uniform_lattice_trajectories_N%d_Con%d_K%d_a0_%.2f_hill%.2f_noise_%.3f_tmax%d_nsim%d_cross1', ...
        N, Con, K, a0, hill, noise, tmax, nsim), '.', 'p');
    load(fullfile(path, fname_in));
    %frac_final(i) = sum(x_out_final > fp(2))/nsim;
    frac_cross1(i) = mean(x_cross1);%sum(x_out_max > fp(2))/nsim;
    
    %
    fname_in = strrep(sprintf('uniform_lattice_trajectories_N%d_Con%d_K%d_a0_%.2f_hill%.2f_noise_%.3f_tmax%d_nsim%d_cross2', ...
        N, Con, K, a0, hill, noise, tmax, nsim), '.', 'p');
    load(fullfile(path, fname_in));
    %frac_final(i) = sum(x_out_final > fp(2))/nsim;
    frac_cross2(i) = mean(x_cross2);%sum(x_out_max > fp(2))/nsim;
    %}
end
%%
h5 = figure(5);
semilogx(noise_all, frac_cross1, 'o--', 'LineWidth', 2, 'MarkerSize', 8);
hold on;
semilogx(noise_all, frac_cross2, 'o--', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('noise level \alpha');
ylabel('crossing probability');
set(gca, 'FontSize', 24);

% figures
fname_str = strrep(sprintf('unif_lattice_crossing_prob_N%d_Con%d_K%d_a0_%.2f_hill%.2f_noise_%.3f_tmax%d_nsim%d', ...
    N, Con, K, a0, hill, noise, tmax, nsim), '.', 'p');
save_fig = 1;
if save_fig > 0
    out_file = fullfile(path, fname_str);
    save_figure_pdf(h5, 10, 8, out_file);
    save_figure_eps(h5, 10, 8, out_file);
end
    
% trick to set log scale
%ax = get(gca);
%xlim(log(ax.XLim));
%set(gca,'XTick',log(ax.XTick),'XTickLabel',ax.XTickLabel);
%}