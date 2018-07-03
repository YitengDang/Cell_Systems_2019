function [ptsum, pt, pEn] = transition_prob_two_signals(n, N, M_int, a0, fN, gN, Rcell, K, Con)
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
% Details of the filename
n_dig = 4; % number of decimal digits used in the filename
mult_dig = 10^n_dig;

% String of the parameters
a0_s = sprintf('%.4d', mult_dig*round(a0, n_dig));
R_s = sprintf('%.4d', mult_dig*round(Rcell, n_dig));
K_s = sprintf('%.4d', mult_dig*round(K, n_dig));
Son_s = sprintf('%.4d', mult_dig*round(Con, n_dig));
%I_s = sprintf('%.4d', mult_dig*round(I, n_dig));
%fname_str = sprintf('tmat_N%d_n%d_a0_%s_R_%s_K_%s_Son_%s_M_int_%d.mat', ...
%    N, n, a0_s, R_s, K_s, Son_s, M_int);
fname_str = 'temp';

%folder = fullfile('..', '..', 'data', 'transition_matrix', 'tmat_n');
folder = 'H:\My Documents\Multicellular automaton\data\two_signals\transition_matrix\tmat_n';
fname = fullfile(folder, fname_str);
if exist(fname,'file') == 2
    tmp = load(fname);
    ptsum = tmp.ptsum;
    pt = tmp.pt;
    pEn = tmp.pEn;
else
    % initial number/fraction of ON cells
    
    p = n/N;
    
    muON = Con + fN.*(p.*Con+(1-p));
    muOFF = 1 + fN.*(p.*Con+(1-p));
    sigmap = sqrt(p.*(1-p).*gN).*(Con-1);

    % vectorized
    muON_mat = repmat(muON, 2, 1);
    sigmap_mat = repmat(sigmap, 2, 1);
    ponon = prod(M_int.*normcdf(M_int.*(-K + muON_mat)./sigmap_mat) + (1-abs(M_int)), 2)';

    muOFF_mat = repmat(muOFF, 2, 1);
    poffoff = prod(M_int.*normcdf(M_int.*(-K + muOFF_mat)./sigmap_mat) + (1-abs(M_int)), 2)';

    %pEn = ponon.^(n').*poffoff.^(N-n');

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