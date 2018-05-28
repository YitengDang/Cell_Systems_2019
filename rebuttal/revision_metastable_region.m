%% Attempts at finding conditions for metastability
clear all
close all
clc
%% Parameters
gridsize = 11;
N=gridsize^2;
a0 = 0.5;
Rcell = 0.2*a0;
K = 16;
Con = 8;

% load fN
[dist, pos] = init_dist_hex(gridsize, gridsize);
dist_vec = dist(1,:);
r = a0*dist_vec(dist_vec>0); % exclude self influence
fN = sum(Rcell*sum(exp(Rcell-r)./r)); % calculate signaling strengthN
%gN = sum(sum((sinh(Rcell)*exp(Rcell-r)./r).^2)); % calculate signaling strength

% functions
fa0 = Rcell*sum(exp(Rcell-a0)./a0);
%mm = @(p) (6 - fN/f(a0))*p + (K-1-fN)/((Con-1)*f(a0));
%mp = @(p) (6 - fN/f(a0))*p + (K-Con-fN)/((Con-1)*f(a0));
%I1 = @(p) p./(1-p) + (K-Con-fN)/((1-p)*(Con-1)*fN);
%I2 = @(p) (1 - 6*f(a0)/fN)*p - (K-1-fN)/((Con-1)*fN);

%% Plot stable region in p, mi space
% values to plot
pv = (0:N)/N;
mOFF = (6 - fN/fa0)*pv + (K-1-fN)/((Con-1)*fa0);
mON = (6 - fN/fa0)*pv + (K-Con-fN)/((Con-1)*fa0);

% truncate values
mOFF(mOFF < 0) = 0;
mOFF(mOFF > 6) = 6;
mON(mON < 0) = 0;
mON(mON > 6) = 6;
mON = floor(mON);
mOFF = floor(mOFF);

% figure
h=figure();
hold on

% transition points
mOFF_trans = zeros(6, 1);
mON_trans = zeros(6, 1);
for i=1:6
    if ~isempty( find(mOFF == 6-i, 1) )
        mOFF_trans(i, 1) = find(mOFF == 6-i, 1) - 1;
    end
    if ~isempty( find(mON == 6-i, 1) )
        mON_trans(i, 1) = find(mON == 6-i, 1) - 1;
    end
    %plot([mOFF_trans(i) mOFF_trans(i)], [0 6-i+1], 'b--', 'LineWidth', 1.2);
    %plot([mON_trans(i) mON_trans(i)], [0 6-i+1], 'r--', 'LineWidth', 1.2);
    plot([mOFF_trans(i)/N mOFF_trans(i)/N], [0 6-i+1], 'b--', 'LineWidth', 1.2);
    plot([mON_trans(i)/N mON_trans(i)/N], [0 6-i+1], 'r--', 'LineWidth', 1.2);
end
plot( (0:N)/N, mON, 'ro', 'LineWidth', 1.2);
plot( (0:N)/N, mOFF, 'bo', 'LineWidth', 1.2);
%xlabel('N_{ON}');
xlabel('p');
ylabel('m_i');
set(gca, 'FontSize', 24);
xlim([0 1]);

qsave = 1;
if qsave 
    fname_str = strrep(sprintf('N%d_a0_%.1f_K_%d_Con_%d', ...
        N, a0, K, Con), '.', 'p');
    fname = fullfile(pwd, 'rebuttal', 'pattern_classification', strcat('nearest_nb_phase_diagram_vs_p_', fname_str));
    save_figure_pdf(h, 10, 8, fname);
    close all
    save(fname); % also save data
end
%% Plot limits for I
I1 = pv./(1-pv) + (K-Con-fN)./((1-pv)*(Con-1)*fN);
I2 = (1 - 6*fa0/fN)*pv - (K-1-fN)/((Con-1)*fN);

figure();
hold on
plot(pv, I1, 'b-');
plot(pv, I2, 'r-');
xlim([0 1]);
ylim([-1 1]);