%% Calculates the mutual information between pin and pout for given pin-pout map
clear all
close all
clc
set(0,'defaulttextinterpreter', 'latex');
%% Parameters
gridsize = 15;
N=gridsize^2;
%Con = 16;
%K = 8;
a0 = 1.5;
Rcell = 0.2*a0;
Conlist = 1:30;
Klist = 1:30;
I_all = zeros(numel(Conlist), numel(Klist));

%calculate fN
[pos,ex,ey] = init_cellpos_hex(gridsize,gridsize);
dist = dist_mat(pos,gridsize,gridsize,ex,ey);
dist_vec = a0*dist(1,:);
r = dist_vec(dist_vec>0); % exclude self influence
fN = sum(sinh(Rcell)*sum(exp(Rcell-r)./r)); % calculate signaling strength

for i=1:numel(Conlist)
    Con = Conlist(i);
    for j=1:numel(Klist)
        K = Klist(j);
        fprintf('Con = %d, K = %d \n', Con, K);
        % Load pin-pout data
        fname_str = sprintf('pin_pout_N%d_Con_%d_K_%d_a0_%d', N, Con, K, a0*10);
        fname = fullfile(pwd, 'data', 'pin_pout', 'analytical', fname_str);
        load(fname, 'prob');
 
        %% Compute mutual information
        % Extreme cases
        %prob = diag(ones(N+1,1)); % autonomy
        %prob = zeros(N+1); prob(1,:)=1; % all OFF 

        % Compute S[ P(p_out) ]
        prob_pout = sum(prob, 2)/(N+1); % P(p_out)
        idx = prob_pout>0;
        S_pout = - sum(prob_pout(idx).*log2(prob_pout(idx)));

        % Compute S[ P(p_out | p_in) ]
        P_pin = 1/(N+1);
        S_pout_pin = zeros(N+1, 1);
        for idx=1:N
            P_pout_pin = prob(:, idx);
            idx2 = P_pout_pin > 0;
            S_pout_pin(idx) = - sum(P_pout_pin(idx2).*log2(P_pout_pin(idx2)));
        end

        % Compute I(pin, pout)
        I = S_pout -sum(P_pin.*S_pout_pin);
        I_scaled = I/log2(N+1); % I rescaled to lie between 0 and 1
        I_all(i,j) = I;
        %figure();
        %plot(prob_pout);

        %% Plot pin-pout map with I
        %{
        set(0, 'defaulttextinterpreter', 'latex');
        p = (0:N)/N;
        h1=figure(1);
        im_fig = imagesc(p, p, prob);
        c = colorbar;
        set(gca, 'YDir', 'Normal','FontSize', 24)
        xlabel('$$p_{in}$$');
        ylabel('$$p_{out}$$');
        title(sprintf('$$I=%.3f = %.3f I_{max}$$', I, I_scaled));
        ylabel(c, '$$P(p_{out}|p_{in})$$', 'Interpreter', 'latex');
        % set invisible parts where count is zero
        set(im_fig, 'AlphaData', prob > 0);
        
        % save figure
        qsave = 0;
        if qsave
            fname_str = sprintf('pin-pout_withI_N%d_a0%d_K_%d_Con_%d',...
                N, 10*a0, K, Con);
            fname = fullfile(pwd, 'figures', 'information', 'all_pin-pout', fname_str); %filename
            save_figure_pdf(h1, 10, 8, fname);
            save_figure_eps(h1, 10, 8, fname);
        end
        %}
   end
end

%% Save mutual information data
qsave = 0;
if qsave
    fname_str = strrep(sprintf('Information_pin-pout_N%d_a0_%.1f_K_%dto%d_Con_%dto%d', ...
        N, a0, Klist(1), Klist(end), Conlist(1), Conlist(end)), '.', 'p');
    fname = fullfile(pwd, 'data', 'pin_pout', 'analytical', strcat(fname_str, '.mat'));
	save(fname, 'I_all', 'Conlist', 'Klist', 'N', 'a0'); 
end
%%
h2 = figure(2);
imagesc(Klist, Conlist, I_all);
xlabel('$$K$$');
ylabel('$$C_{ON}$$');
set(gca, 'YDir', 'Normal','FontSize', 24);
c = colorbar;
ylabel(c, '$$I(p_{out}, p_{in})$$', 'Interpreter', 'latex');
ticks = [1 5 10 15 20 25 30];
set(gca, 'xtick', ticks, 'ytick', ticks);
hold on

%plot([1 30], [4.5 4.5]);
plot([fN+1 fN+1], [Conlist(1) Conlist(end)+1], 'g', 'LineWidth', 1.5) % ON region
plot([Klist Klist(end)+1], [Klist./(1+fN) (Klist(end)+1)/(1+fN)], 'r', 'LineWidth', 1.5) % OFF region
if (Klist(end/2)-fN) < (Klist(end/2)-1)/fN && (Klist(end/2)-1)/fN > Conlist(end/2)
    plot([Klist Klist(end)+1], [Klist Klist(end)+1]-fN, 'k', 'LineWidth', 1.5) % autonomous 1
    plot([Klist Klist(end)+1], ([Klist Klist(end)+1]-1)/fN, 'k', 'LineWidth', 1.5) % autonomous 2
end
%plot(Klist, 2*Klist/(1+fN)-1, '--k', 'LineWidth', 1.5) % 0.5(Con_1)(1+fN) = K

qsave = 1;
if qsave
    fname_str = strrep(sprintf('Information_K_Con_map_N%d_a0_%.1f_infinite_Hill_phases',...
        N, a0), '.', 'p');
    fname = fullfile(pwd, 'figures', 'information', fname_str); %filename
    save_figure_pdf(h2, 10, 8, fname);
    save_figure_eps(h2, 10, 8, fname);
end