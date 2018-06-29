function [ptsum, pt, pEn] = transition_prob_with_I(dist, a0, Rcell, K, Son, n, I)
% Calculate the estimated transition probability matrix for different
% initial number n of ON cells. pEn is the equilibrium probability.

if nargin < 7
    I = 0;
end

dist_vec = dist(1,:);
N = numel(dist_vec);
r = a0*dist_vec(dist_vec>0); % exclude self influence
fN = sum(sinh(Rcell)*sum(exp(Rcell-r)./r)); % calculate signaling strength
gN = sum(sum((sinh(Rcell)^2)*exp(2*(Rcell-r))./(r.^2))); % calculate noise variance strength

% Details of the filename
n_dig = 4; % number of decimal digits used in the filename
mult_dig = 10^n_dig;
% String of the parameters
a0_s = sprintf('%.4d', mult_dig*round(a0, n_dig));
R_s = sprintf('%.4d', mult_dig*round(Rcell, n_dig));
K_s = sprintf('%.4d', mult_dig*round(K, n_dig));
Son_s = sprintf('%.4d', mult_dig*round(Son, n_dig));
I_s = sprintf('%.4d', mult_dig*round(I, n_dig));

fname_str = sprintf('tmat_N%d_n%d_a0_%s_R_%s_K_%s_Con_%s_I_%s.mat', ...
    N, n, a0_s, R_s, K_s, Son_s, I_s);
%folder = fullfile('..', '..', 'data', 'transition_matrix', 'tmat_n');
folder = 'H:\My Documents\Multicellular automaton\data\main\transition_matrix\tmat_n';
fname = fullfile(folder, fname_str);

if exist(fname,'file') == 2
    tmp = load(fname);
    ptsum = tmp.ptsum;
    pt = tmp.pt;
    pEn = tmp.pEn;
else
    
    p = n/N;
    muON = Son + fN*(p*Son+(1-p)+(Son-1)*(1-p)*I);
    muOFF = 1 + fN*(p*Son+(1-p)-(Son-1)*p*I);
    sigmap = sqrt(p.*(1-p)*gN)*(Son-1);
    
    ponon = normcdf((-K+muON)./sigmap);
    poffoff = normcdf((K-muOFF)./sigmap);
    
    yminv = 0:n;
    yplusv = 0:(N-n);
    
    pt = zeros(numel(yminv), numel(yplusv));
    ptsum = zeros(N+1,1);
    for ymin = yminv
        for yplus = yplusv
            tmp = ponon^(n-ymin)*(1-ponon)^ymin*poffoff^(N-n-yplus)*(1-poffoff)^yplus;
            pt(ymin+1,yplus+1) = tmp*nchoosek(n,ymin)*nchoosek(N-n,yplus);
            idx = n + yplus - ymin;
            ptsum(idx + 1) = ptsum(idx + 1) + pt(ymin+1,yplus+1);
        end
    end
    
    pEn = pt(1, 1);
    save(fname, 'ptsum', 'pt', 'pEn');
end