function [ptsum, pt, pEn] = transition_prob(dist, a0, Rcell, K, Son, n)
% Calculate the estimated transition probability matrix for different
% initial number n of ON cells. pEn is the equilibrium probability.

dist_vec = dist(1,:);
N = numel(dist_vec);
r = a0*dist_vec(dist_vec>0); % exclude self influence
fN = sum(sinh(Rcell)*sum(exp(Rcell-r)./r)); % calculate signaling strength
gN = sum(sum((sinh(Rcell)^2)*exp(2*(Rcell-r))./(r.^2))); % calculate noise variance strength

p = n/N;
muON = Son + fN*(p*Son+(1-p));
muOFF = 1 + fN*(p*Son+(1-p));
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