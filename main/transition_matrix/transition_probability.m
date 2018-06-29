% Calculate the estimated transition probability of a certain number of 
% cells activating or deactivating for a specified initial number of ON cells

clear variables
close all

% Parameters of the system
gridsize = 11;
N = gridsize^2;
a0 = 0.5;
Rcell = 0.2*a0;

Son = 8;
K = 16;

% use hexagonal lattice
[pos,ex,ey] = init_cellpos_hex(gridsize,gridsize);
dist = dist_mat(pos,gridsize,gridsize,ex,ey);

dist_vec = dist(1,:);
r = a0*dist_vec(dist_vec>0); % exclude self influence
fN = sum(sinh(Rcell)*sum(exp(Rcell-r)./r)); % calculate signaling strength
gN = sum(sum((sinh(Rcell)^2)*exp(2*(Rcell-r))./(r.^2))); % calculate noise variance strength

%n = 0:N;
n=60;
p = n/N;
muON = Son + fN*(p*Son+(1-p));
muOFF = 1 + fN*(p*Son+(1-p));
sigmap = sqrt(p.*(1-p)*gN)*(Son-1);

ponon = normcdf((-K+muON)./sigmap);
poffoff = normcdf((K-muOFF)./sigmap);

pEn = ponon.^(n).*poffoff.^(N-n);

yminv = 0:n;
yplusv = 0:(N-n);

pt = zeros(numel(yminv), numel(yplusv));
ptn = zeros(N+1,1);
for ymin = yminv
    for yplus = yplusv
        tmp = ponon^(n-ymin)*(1-ponon)^ymin*poffoff^(N-n-yplus)*(1-poffoff)^yplus;
        pt(ymin+1,yplus+1) = tmp*nchoosek(n,ymin)*nchoosek(N-n,yplus);
        idx = n + yplus - ymin;
        ptn(idx + 1) = ptn(idx + 1) + pt(ymin+1,yplus+1);
    end
end
       
imagesc(yminv, yplusv, log(pt'))
c = colorbar;
set(gca,'ydir','normal')
