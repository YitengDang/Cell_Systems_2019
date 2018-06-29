function [ptsum, pt, pEn] = transition_prob(dist, a0, Rcell, K, Son, n, M_int)
% Calculate the estimated transition probability of a certain number of 
% cells activating or deactivating for a specified initial number of ON cells

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
%I_s = sprintf('%.4d', mult_dig*round(I, n_dig));
fname_str = sprintf('tmat_N%d_n%d_a0_%s_R_%s_K_%s_Son_%s_M_int_%d.mat', ...
    N, n, a0_s, R_s, K_s, Son_s, M_int);

%folder = fullfile('..', '..', 'data', 'transition_matrix', 'tmat_n');
folder = 'H:\My Documents\Multicellular automaton\data\main\transition_matrix\tmat_n';
fname = fullfile(folder, fname_str);
if exist(fname,'file') == 2
    tmp = load(fname);
    ptsum = tmp.ptsum;
    pt = tmp.pt;
    pEn = tmp.pEn;
else
    % initial number/fraction of ON cells
    %n=60; 
    p = n/N;

    muON = Son + fN*(p*Son+(1-p));
    muOFF = 1 + fN*(p*Son+(1-p));
    sigmap = sqrt(p.*(1-p)*gN)*(Son-1);
    
    zon = (K-muON)./sigmap;
    zoff = (K-muOFF)./sigmap;
    if M_int==1
        ponon = 1-normcdf(zon);
        poffoff = normcdf(zoff);
    elseif M_int ==-1
        ponon = normcdf(zon);
        poffoff = 1-normcdf(zoff);
    end

    %pEn = ponon.^(n).*poffoff.^(N-n); %already follows from calculation
    %below

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

end