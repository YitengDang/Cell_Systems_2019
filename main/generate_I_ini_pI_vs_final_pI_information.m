% This script uses the calculated P(p_out, I_out | p_in, I_in) to calculate
% mutual information I(p_out, I_out | p_in, I_in)
close all
clear variables
%warning off

% Parameters of the system
gridsize = 15;
N = gridsize^2;
a0 = 0.5;
Rcell = 0.2*a0;
K = 15;
Con = 8;
%a0all = [0.5 1.5 1.5];
%Kall = [15 6 20];
%Conall = [8 15 15];

p1 = floor(0.1*N)/N;
dp = round(0.05*N)/N;
p_all = p1:dp:(1-p1);
I_all = -0.05:0.05:0.4;

% use hexagonal lattice
[dist, ~] = init_dist_hex(gridsize, gridsize);

% Calculate the signaling strength
dist_vec = dist(1,:);
r = a0*dist_vec(dist_vec>0); % exclude self influence
fN = sum(sinh(Rcell)*sum(exp(Rcell-r)./r)); % calculate signaling strength

%% load the calculated map
%uigetfile(fullfile(pwd, 'data', 'pIin_pIout'));
fname_str = strrep(sprintf('pIin_pIout_N%d_Con_%d_K_%d_a0_%.1f', ...
    N, Con, K, a0), '.', 'p');
v = 1;
fname = fullfile(pwd, 'data', 'pIin_pIout', strcat(fname_str,'-v', int2str(v), '.mat'));
load(fname)
close all;

%% Calculation of mutual information 

% Calculate P(p_out, I_out)
p_edges = 0:0.05:1;
I_edges = -0.2:0.05:0.6;
dI = (I_edges(end)-I_edges(1))/numel(I_edges);

[count] = histcounts2(p_out(:), I_out(:), p_edges, I_edges);
prob_out = transpose(count)/sum(sum(count)); 

% Calculate S[P(p_out, I_out)]
idx = prob_out > 0;
S_out = -sum(sum(prob_out(idx).*log(prob_out(idx))));

% Calculate S[P(p_out, I_out | p_in, I_in)]
P_in = 1/(numel(p_all)*numel(I_all));
S_out_in_temp = zeros(numel(p_all), numel(I_all));
for i=1:numel(p_all)
    for j=1:numel(I_all)
        [count2] = histcounts2(p_out(i,j, :), I_out(i,j,:), p_edges, I_edges);
        prob_out_in = transpose(count2)/sum(sum(count2));
        idx2 = prob_out_in > 0;
        S_out_in_temp(i, j) = - sum(prob_out_in(idx2).*log2(prob_out_in(idx2)));
    end
end
S_out_in = sum(sum(P_in.*S_out_in_temp));
I_out_in = S_out - S_out_in;
%}

%% Plot P(p_out, I_out)
set(0, 'defaulttextinterpreter', 'latex');
h1=figure(1);
hold on
im_fig=imagesc(p_edges, I_edges, prob_out);
set(im_fig, 'AlphaData', count' > 0);
set(gca, 'YDir', 'Normal','FontSize', 24);
c = colorbar;
xlabel('$$p_{out}$$');
ylabel('$$I_{out}$$');
ylabel(c, '$$P(p_{out}, I_{out})$$', 'interpreter', 'latex');
title(sprintf('$$I(p_{out}, I_{out}|p_{in}, I_{in}) \\approx %.4f$$ bits', I_out_in));
xlim([0-dp/2 1+dp/2]);
ylim([I_edges(1)-dI I_edges(end)+dI]);

% Plot initial p, I as grid 
% values are at center of each rectangle of grid
[X,Y] = meshgrid([p_all-dp/2 p_all(end)+dp/2], I_all(1)-dp/2:0.01:I_all(end)+dp/2 );
plot(X,Y,'r.');

[X,Y] = meshgrid(p_all(1)-dp/2:0.01:p_all(end)+dp/2, [I_all-dI/2 I_all(end)+dI/2]);
plot(X,Y,'r.');

%Imax = log2();
%fprintf('I = %.3f Imax', ); 

% save figure
qsave = 0;
if qsave
    fname_str = strrep(sprintf('pI_out_withI_N%d_a0%d_K_%d_Con_%d_pout_%.3fto%.3f_Iout_%.3fto%.3f',...
        N, 10*a0, K, Con, p_edges(1), p_edges(end), I_edges(1), I_edges(end)), '.', 'p');
    fname = fullfile(pwd, 'figures', 'pIin_pIout', fname_str); %filename
    save_figure_pdf(h1, 10, 8, fname);
    save_figure_eps(h1, 10, 8, fname);
end