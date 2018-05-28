% Visualize the dynamics of a saved data

clear variables
close all
warning off

% Times to analyse
t_vec = [0 10];

% Save figure?
save_fig = 1;

data_path = 'D:\eduardopavinat\Dropbox\Matlab codes\data_onecelltype_entropy';
[fname, path, ~] = uigetfile(fullfile(data_path,'dynamics_noise'));
load(fullfile(path,fname))
[~,name,~] = fileparts(fname);

% First neighbor distance
eps = 1e-5;
dist_vec = get_unique_distances(dist, eps);
dist1 = dist_vec(2);
idx = 1*(dist < dist1+eps & dist > dist1-eps);

for i = 1:numel(t_vec)
    hin = figure(1);
    update_cell_figure(hin, pos, a0, cells_hist{t_vec(i)+1}, cell_type, t_vec(i))
    if save_fig > 0
        out_file = fullfile(pwd, 'figures', 'dynamics_nonoise', ...
            strcat(name,'_t_',int2str(t_vec(i))));
        save_figure_pdf(hin, 10, 10, out_file);
    end
    hhist = figure(4);
    nn = idx*cells_hist{t_vec(i)+1};
    subplot(2,1,1)
    histogram(nn(cells_hist{t_vec(i)+1} == 1), 0:7)
    set(gca, 'FontSize', 20, 'XTick', 0:6)
    title(sprintf('p = %.2f, t = %d', sum(cells_hist{t_vec(i)+1})/N, t_vec(i)), ...
        'FontSize', 20)
    ylabel('Count', 'FontSize', 24)
    xlabel('m_{i=ON}', 'FontSize', 24)
    xlim([0 7])
    subplot(2,1,2)
    histogram(nn(cells_hist{t_vec(i)+1} == 0), 0:7)
    set(gca, 'FontSize', 20, 'XTick', 0:6)
    ylabel('Count', 'FontSize', 24)
    xlabel('m_{i=OFF}', 'FontSize', 24)
    xlim([0 7])
    if save_fig > 0
        out_file = fullfile(pwd, 'figures', 'dynamics_nonoise', ...
            strcat(name,'_hist_t_',int2str(t_vec(i))));
        save_figure_pdf(hhist, 7, 8, out_file);
    end    
end

fa0 = sinh(Rcell)*sum(exp(Rcell-a0)./a0);
for i = 1:numel(cells_hist)
    nn = idx*cells_hist{i};
    p = sum(cells_hist{i})/N;
    p_vec(i) = p;
    nn_on = mean(nn(cells_hist{i} == 1));
    nn_off = mean(nn(cells_hist{i} == 0));
    if p == 1 || p==0
        I_test(i) = 0;
    else
        I_test(i) = fa0/fN*max((nn_on - 6*p)/(1-p), (6*p - nn_off)/p);
    end
end
h9 = figure(9);
hold on
plot(0:t, I_test, 'b-o', 'Linewidth', 1.5)
plot(0:t, I, 'r-o', 'Linewidth', 1.5)
hold off
legend({'Approx.', 'Exact'}, 'Location', 'northwest')
set(gca, 'FontSize', 20)
ylabel('Spatial order I', 'FontSize', 24)
xlabel('Time (steps)', 'FontSize', 24)
box on
if save_fig > 0
    out_file = fullfile(pwd, 'figures', 'dynamics_nonoise', strcat(name,'_Itest'));
    save_figure_pdf(h9, 7, 6, out_file);
end

h2 = figure(2);
subplot(2,1,1)
plot(0:t, I, 'r-o', 'Linewidth', 1.5)
xlim([0 t])
set(gca, 'FontSize', 20, 'XTickLabel', '')
%xlabel('Time (steps)', 'FontSize', 24)
ylabel('Moran index', 'FontSize', 24)
% if save_fig > 0
%     out_file = fullfile(pwd, 'figures', 'dynamics_nonoise', strcat(name,'_I'));
%     save_figure_pdf(h2, 8, 6, out_file);
% end

subplot(2,1,2)
plot(0:t, Non/N, 'r-o', 'Linewidth', 1.5)
xlim([0 t])
set(gca, 'FontSize', 20)
xlabel('Time (steps)', 'FontSize', 24)
ylabel('p = N_{ON}/N', 'FontSize', 24)
if save_fig > 0
    out_file = fullfile(pwd, 'figures', 'dynamics_nonoise', strcat(name,'_p'));
    save_figure_pdf(h2, 8, 10, out_file);
end

h3 = figure(3);
plot(0:t, -mom/N, 'r-o', 'Linewidth', 1.5)
set(gca, 'FontSize', 20)
xlabel('Time (steps)', 'FontSize', 24)
ylabel('H/N', 'FontSize', 24)
xlim([0 t])

if save_fig > 0
    out_file = fullfile(pwd, 'figures', 'dynamics_nonoise', strcat(name,'_H'));
    save_figure_pdf(h3, 8, 6, out_file);
end

h6 = figure(6);
hold on
piv = (0:N)/N;
Iv = -0.2:0.1:1;
E = @(p, I) -0.5*(Son-1)*(1 + 4*fN.*p.*(1-p).*I + fN*(2*p-1).^2) ...
    -(2*p-1).*(0.5*(Son+1)*(1+fN) - K);
u = @(p, I) 2*(Son-1)*fN*(2*p-1).*(1-I)+(Son+1)*(1+fN)-2*K;
v = @(p, I) 2*(Son-1)*fN*p.*(1-p);
[p_i,Imesh] = meshgrid(piv, Iv);
contourf(piv, Iv, E(p_i, Imesh), 'LineStyle', 'none')
colormap('summer')
[p_i,Im] = meshgrid(piv(1:4:end), Iv);

fa0 = sinh(Rcell)*exp(Rcell-a0)/a0;
Imax = @(p) (6 - 4./sqrt(p*N) - 6./(p*N) - 6*p)*fa0./(1-p)/fN;
Imax2 = @(p) (6*p - 4./sqrt((1-p)*N) - 6./((1-p)*N))*fa0./p/fN;
%plot(piv, Imax(piv), 'b')
plot(piv, max(Imax(piv), Imax2(piv)), '--b', 'Linewidth', 1.5)
ylim([-0.05 1])

plot(p_vec, I, 'k-o')
hold off
