% This script uses the calculated P(p_out, I_out | p_in, I_in) to perform
% Bayesian inference and calculate P(p_in, I_in | p_out, I_out) maps
close all
clear variables
%warning off

% Parameters of the system
gridsize = 15;
N = gridsize^2;
%a0 = 0.5;
%Rcell = 0.2*a0;
%K = 15;
%Con = 8;
a0all = [0.5 1.5 1.5];
Kall = [15 6 20];
Conall = [8 15 15];
var_idx = 3;
a0 = a0all(var_idx);
K = Kall(var_idx);
Con = Conall(var_idx);
Rcell = 0.2*a0;

p1 = floor(0.1*N)/N;
dp = round(0.05*N)/N;
p_all = p1:dp:(1-p1);
I_all = -0.05:0.05:0.4;

% use hexagonal lattice
[dist, pos] = init_dist_hex(gridsize, gridsize);

% filename of the saved data
fname_str = strrep(sprintf('pIin_pIout_N%d_Con_%d_K_%d_a0_%.1f', ...
    N, Con, K, a0), '.', 'p');

% Calculate the signaling strength
dist_vec = dist(1,:);
r = a0*dist_vec(dist_vec>0); % exclude self influence
fN = sum(sinh(Rcell)*sum(exp(Rcell-r)./r)); % calculate signaling strength

%% load the calculated map
%uigetfile(fullfile(pwd, 'data', 'pIin_pIout'));
v = 1;
fname = fullfile(pwd, 'data', 'pIin_pIout', strcat(fname_str,'-v', int2str(v), '.mat'));
load(fname)
close all;
%% Plot P(p_out, I_out | p_in, I_in) map
I_in = 0.2;
%p_in = 0.5;
for idx = 1:numel(p_all)
    p_in = p_all(idx);
    
    idx_p = find(round(p_all, 10) == round(round(p_in/dp)*dp, 10), 1);
    idx_I = find(I_all == I_in, 1);

    if isempty(idx_p) || isempty(idx_I)
        disp('Error, enter valid p and I!');
    end

    % histograms with only p or only I
    %figure();
    %histogram(p_out(idx_p, idx_I, :));
    %figure();
    %histogram(I_out(idx_p, idx_I, :));

    % heat maps with P(p_out | p_in, I_in) and P(I_out | p_in, I_in)
    pedges = 0:0.05:1;
    Iedges = -0.2:0.05:0.5;
    [counts, ~, ~] = histcounts2(I_out(idx_p, idx_I, :), p_out(idx_p, idx_I, :), Iedges, pedges);
    % Note I: 1st index vert., p: 2nd index horz.

    %h1=figure(1);
    %figure()
    %imagesc(pedges, Iedges, counts);
    %title(sprintf('p = %.3f', p_in));
end

%% plot expected p_out given (p_in, I_in)
% calculate expected final (p,I) for each initial (p,I)
p_out_av = transpose(mean(p_out, 3));
I_out_av = transpose(mean(I_out, 3));
set(0, 'defaulttextinterpreter', 'latex');

% figure
h1 = figure(1);
imagesc(p_all, I_all, p_out_av); % plot final states
set(gca,'Ydir','normal','FontSize', 24)
% set colorbar and labels
c = colorbar;
set(gca, 'FontSize', 24);
xlabel('$$p_{in}$$');
ylabel('$$I_{in}$$');
ylabel(c,'$$\langle p_{out} \rangle$$', 'interpreter', 'latex');

% save
qsave=0;
if qsave
   fname_str = strrep(sprintf('pIin_mean_p_out_N%d_a0_%.1f_Con%d_K%d_I_%.3f_%.3f_p_%.3f_%.3f',...
       N, a0, Con, K, I_all(1), I_all(end), p_all(1), p_all(end)), '.', 'p');
   fname = fullfile(pwd, 'figures', 'pIin_pIout', fname_str);
   save_figure_pdf(h1, 10, 8, fname);
   save_figure_eps(h1, 10, 8, fname);
end

%% plot expected I_out given (p_in, I_in)
h2 = figure(2);
imagesc(p_all, I_all, I_out_av); % plot final states
set(gca,'Ydir','normal','FontSize', 24)
% set colorbar and labels
c = colorbar;
set(gca, 'FontSize', 24);
xlabel('$$p_{in}$$');
ylabel('$$I_{in}$$');
ylabel(c,'$$\langle I_{out} \rangle$$', 'interpreter', 'latex');

% save
qsave=0;
if qsave
   fname_str = strrep(sprintf('pIin_mean_I_out_N%d_a0_%.1f_Con%d_K%d_I_%.3f_%.3f_p_%.3f_%.3f',...
       N, a0, Con, K, I_all(1), I_all(end), p_all(1), p_all(end)), '.', 'p');
   fname = fullfile(pwd, 'figures', 'pIin_pIout', fname_str);
   save_figure_pdf(h2, 10, 8, fname);
   save_figure_eps(h2, 10, 8, fname);
end


%% Plot the average number of steps it takes to reach equilibrium
h3 = figure(3);
%tav = zeros(numel(p_all), numel(I_all));
imagesc(p_all, I_all, t_av');
set(gca,'Ydir','normal','FontSize', 24)
c = colorbar;
xlabel('p', 'FontSize', 24)
ylabel('I', 'FontSize', 24)
ylabel(c, '$$ \langle t_{eq} \rangle$$', 'FontSize', 24, 'Interpreter', 'latex')

qsave=0;
if qsave
   fname_str = strrep(sprintf('pIin_pIout_N%d_a0_%.1f_Con%d_K%d_I_%.3f_%.3f_p_%.3f_%.3f_tav',...
       N, a0, Con, K, I_all(1), I_all(end), p_all(1), p_all(end)), '.', 'p');
   fname = fullfile(pwd, 'figures', 'pIin_pIout', fname_str);
   save_figure_pdf(h3, 10, 8, fname);
   save_figure_eps(h3, 10, 8, fname);
end

%% Plot P(p_out, I_out)
% Calculate P(p_out, I_out)
p_edges = 0:0.05:1;
I_edges = -0.2:0.05:0.6;
dI = (I_edges(end)-I_edges(1))/numel(I_edges);

[count] = histcounts2(p_out(:), I_out(:), p_edges, I_edges);
prob_out = transpose(count)/sum(sum(count)); 
%% Plot P(P_out, I_out)
h4 = figure(4);
hold on
% Plot P(p_out, I_out)
set(0, 'defaulttextinterpreter', 'latex');
im_fig=imagesc(p_edges, I_edges, prob_out);
set(im_fig, 'AlphaData', count' > 0);
set(gca, 'YDir', 'Normal','FontSize', 24);
c = colorbar;
xlabel('$$p_{out}$$');
ylabel('$$I_{out}$$');
ylabel(c, '$$P(p_{out}, I_{out})$$', 'interpreter', 'latex');
%title(sprintf('$$I(p_{out}, I_{out}|p_{in}, I_{in}) \\approx %.4f$$ bits', I_out_in));
xlim([0-dp/2 1+dp/2]);
ylim([I_edges(1)-dI I_edges(end)+dI]);

% Plot initial p, I as grid 
% values are at center of each rectangle of grid
[X,Y] = meshgrid([p_all-dp/2 p_all(end)+dp/2], I_all(1)-dp/2:0.01:I_all(end)+dp/2 );
plot(X,Y,'r.');

[X,Y] = meshgrid(p_all(1)-dp/2:0.01:p_all(end)+dp/2, [I_all-dI/2 I_all(end)+dI/2]);
plot(X,Y,'r.');

qsave = 0;
if qsave
   fname_str = strrep(sprintf('pIin_pIout_N%d_a0_%.1f_Con%d_K%d_I_%.3f_%.3f_p_%.3f_%.3f_out_map',...
       N, a0, Con, K, I_all(1), I_all(end), p_all(1), p_all(end)), '.', 'p');
   fname = fullfile(pwd, 'figures', 'pIin_pIout', fname_str);
   save_figure_pdf(h4, 10, 8, fname);
   save_figure_eps(h4, 10, 8, fname);
end
%% Plot P(p_out | p_in)
p_count = zeros(numel(p_all), numel(p_edges)-1);

for idx_p=1:numel(p_all)
%idx_p = 1;
p_count(idx_p,:) = histcounts( p_out(idx_p, :, :), pedges);
end

P_pout_pin = transpose( p_count./sum(p_count,2) );
pout_bins = (p_edges(1:end-1)+p_edges(2:end))/2;

h5 = figure(5);
im_fig = imagesc(p_all, pout_bins, P_pout_pin);
set(gca,'Ydir','normal','FontSize', 20)
% set invisible parts where count is zero
set(im_fig, 'AlphaData', p_count' > 0);
% set colorbar and labels
c = colorbar;
%c.Label.String = 'Probability';
xlabel('$$p_{in}$$', 'FontSize', 24)
ylabel('$$p_{out}$$', 'FontSize', 24)
ylim(c, [0 1]);
ylabel(c, '$$P(p_{out}|p_{in})$$', 'FontSize', 24, 'Interpreter', 'latex')

qsave = 1;
if qsave
   fname_str = strrep(sprintf('pin_pout_from_pI_sims_N%d_a0_%.1f_Con%d_K%d_I_%.3f_%.3f_p_%.3f_%.3f_map',...
       N, a0, Con, K, I_all(1), I_all(end), p_all(1), p_all(end)), '.', 'p');
   fname = fullfile(pwd, 'figures', 'pIin_pIout', fname_str);
   save_figure_pdf(h5, 10, 8, fname);
   save_figure_eps(h5, 10, 8, fname);
end

%% Plot P(p_out | p_in)
I_count = zeros(numel(I_all), numel(I_edges)-1);

for idx_I=1:numel(I_all)
    %idx_p = 1;
    I_count(idx_I, :) = histcounts( I_out(:, idx_I, :), I_edges);
end

P_Iout_Iin = transpose( I_count./sum(I_count, 2));
Iout_bins = (I_edges(1:end-1)+I_edges(2:end))/2;

h6 = figure(6);
im_fig = imagesc(I_all, Iout_bins, P_Iout_Iin);
%im_fig = imagesc([0 1], [1 2], P_Iout_Iin);
set(gca,'Ydir','normal','FontSize', 20)
% set invisible parts where count is zero
set(im_fig, 'AlphaData', I_count' > 0);
% set colorbar and labels
c = colorbar;
%c.Label.String = 'Probability';
xlabel('$$I_{in}$$', 'FontSize', 24)
ylabel('$$I_{out}$$', 'FontSize', 24)
ylabel(c, '$$P(I_{out}|I_{in})$$', 'FontSize', 24, 'Interpreter', 'latex')
ylim(c, [0 1]);

qsave = 1;
if qsave
   fname_str = strrep(sprintf('Iin_Iout_from_pI_sims_N%d_a0_%.1f_Con%d_K%d_I_%.3f_%.3f_p_%.3f_%.3f_map',...
       N, a0, Con, K, I_all(1), I_all(end), p_all(1), p_all(end)), '.', 'p');
   fname = fullfile(pwd, 'figures', 'pIin_pIout', fname_str);
   save_figure_pdf(h6, 10, 8, fname);
   save_figure_eps(h6, 10, 8, fname);
end