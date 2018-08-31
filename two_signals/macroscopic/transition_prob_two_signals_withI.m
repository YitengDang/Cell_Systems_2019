function [ptsum, pt, pEn] = transition_prob_two_signals_withI(n, I, N, M_int, a0, fN, gN, rcell, K, Con, Coff)
% Calculate the estimated transition probability of a certain number of 
% cells activating or deactivating for a specified initial number of ON cells
%{
dist_vec = dist(1,:);
N = numel(dist_vec);
r = a0*dist_vec(dist_vec>0); % exclude self influence
fN = sum(sinh(Rcell)*sum(exp(Rcell-r)./r)); % calculate signaling strength
gN = sum(sum((sinh(Rcell)^2)*exp(2*(Rcell-r))./(r.^2))); % calculate noise variance strength
%}
%%
% String of the parameters
a0_s = sprintf('%.2f', a0);
R_s = sprintf('%.2f', rcell);
K_s = sprintf('%d_%d_%d_%d', K(1,1), K(1,2), K(2,1), K(2,2));
Con_s = sprintf('%d_%d', Con(1), Con(2));
M_int_s = sprintf('%d_%d_%d_%d', M_int(1,1), M_int(1,2), M_int(2,1), M_int(2,2));
%I_s = sprintf('%d', mult_dig*round(I, n_dig));
fname_str = strrep(sprintf('tmat_N%d_n1_%d_n2_%d_a0_%s_rcell_%s_K_%s_Con_%s_M_int_%s',...
    N, n(1), n(2), a0_s, R_s, K_s, Con_s, M_int_s), '.', 'p');
%fname_str = 'temp';
%folder = fullfile('..', '..', 'data', 'transition_matrix', 'tmat_n');
folder = 'H:\My Documents\Multicellular automaton\data\two_signals\transition_matrix\tmat_n';
fname = fullfile(folder, strcat(fname_str, '.mat'));
%%
if exist(fname,'file') == 2
    tmp = load(fname);
    ptsum = tmp.ptsum;
    pt = tmp.pt;
    pEn = tmp.pEn;
else
    % initial number/fraction of ON cells
    p = n/N;
    
    Y_nei_avg = fN.*(p.*Con+(1-p) + (Con-1).*(1-p).*I); % neighbour contributions
    mu_mat = repmat(Y_nei_avg, 2, 1);

    sigma = sqrt(p.*(1-p).*gN).*(Con-1);
    sigma_mat = repmat(sigma, 2, 1);

    % Calculate self-contributions
    %off-diagonal (mean-field) self-contributions, other gene is unknown
    self_off_diag = zeros(2); 
    self_off_diag(1,2) = (Con(2)*p(2) + Coff(2)*(1-p(2)));
    self_off_diag(2,1) = (Con(1)*p(1) + Coff(1)*(1-p(1)));

    Y_self_ON_mat = diag(Con) + self_off_diag;
    Y_self_OFF_mat = diag(Coff) + self_off_diag;

    % --- Calculate Ponon, Poffoff ---
    % vectorized
    ponon = prod(1 - abs(M_int).*normcdf(0, M_int.*(mu_mat + Y_self_ON_mat - K), sigma_mat), 2 )';
    poffon = prod(1 - abs(M_int).*normcdf(0, M_int.*(mu_mat + Y_self_OFF_mat - K), sigma_mat), 2)';
    poffoff = 1 - poffon;

    % Calculate pt, ptn
    ptsum = zeros(N+1);
    yminv1 = 0:n(1);
    yminv2 = 0:n(2);
    yplusv1 = 0:(N-n(1));
    yplusv2 = 0:(N-n(2));
    pt = zeros(numel(yminv1), numel(yplusv1), numel(yminv2), numel(yplusv2));
    for ymin1 = yminv1
        for yplus1 = yplusv1
            for ymin2 = yminv2
                for yplus2 = yplusv2
                    %disp([ymin1 yplus1 ymin2 yplus2]);
                    ymin = [ymin1 ymin2];
                    yplus = [yplus1 yplus2];
                    tmp = ponon.^(n-ymin).*(1-ponon).^ymin.*poffoff.^(N-n-yplus).*(1-poffoff).^yplus;
                    tmp(1) = tmp(1)*nchoosek(n(1),ymin1)*nchoosek(N-n(1),yplus1);
                    tmp(2) = tmp(2)*nchoosek(n(2),ymin2)*nchoosek(N-n(2),yplus2);
                    pt(ymin1+1,yplus1+1,ymin2+1,yplus2+1) = prod(tmp);
                    idx = n + yplus - ymin;
                    ptsum(idx(1)+1, idx(2)+1) = ptsum(idx(1)+1, idx(2)+1) + prod(tmp);
                end
            end
        end
    end
    pEn = pt(1, 1, 1, 1);
    save(fname, 'ptsum', 'pt', 'pEn');
end

end