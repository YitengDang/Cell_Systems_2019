function entropy = entropy_two_cell(dist_vec, Nv, Son, K, a0, Rcell)
% Calculate the analytical entropy given the parameters of the system
% dist_vec is a vector with all the distances r_i. dist_vec has N elements.
% omega k gives the number of eq. states for a given fraction of ON cells.
% Nv, Son and K are 2x1 vectors with the parameters for each cell type

dist_vec = dist_vec*a0;
N = sum(Nv);
r = dist_vec(dist_vec>0); % exclude self influence
fN = sum(sinh(Rcell)*sum(exp(Rcell-r)./r)); % calculate signaling strength
gN = sum(sum((sinh(Rcell)^2)*exp(2*(Rcell-r))./(r.^2))); % calculate noise variance strength

% nchoosek calculation is saved in data to reduce computation time. Then it
% is only calculated if the file does not exists yet. It saves the values
% of all Binom(N,n) for a given N and n = 0:N. The values are then used in
% the vector nck.
nck_all = {[],[]};
i = 0;
p = {[], []};
for Nl = Nv
    i = i+1;
    filename = fullfile(pwd,'data','nchoosek',sprintf('%d.mat', Nl));
    if exist(filename, 'file') == 2
        tmp = load(filename);
        nck_all{i} = tmp.nck;
    else
        nck = zeros(Nl+1, 1);
        for n = 0:Nl
            nck(n+1) = nchoosek(Nl,n);
        end
        nck_all{i} = nck;
        save(filename, 'nck');
    end
    p{i} = (0:Nl)/Nl;
end

% Use meshgrid to calculate both fractions in a single shot (avoid slow
% loops in Matlab)
[nck1, nck2] = meshgrid(nck_all{1}, nck_all{2});
[p1, p2] = meshgrid(p{1}, p{2});

mu = fN/N*(Nv(1)*((Son(1)-1)*p1 + 1) + Nv(2)*((Son(2)-1)*p2 + 1));
sigma = sqrt(gN/N*(Nv(1)*(Son(1) - 1)^2*p1.*(1-p1) +  ...
    Nv(2)*(Son(2) - 1)^2*p2.*(1-p2)));
Ponon1 = normcdf((mu+Son(1)-K(1))./sigma);
Ponon2 = normcdf((mu+Son(2)-K(2))./sigma);
Poffoff1 = normcdf((K(1)-mu-1)./sigma);
Poffoff2 = normcdf((K(2)-mu-1)./sigma);
omega = sum(sum(Ponon1.^(Nv(1)*p1).*Ponon2.^(Nv(2)*p2) ...
    .*Poffoff1.^(Nv(1)*(1-p1)).*Poffoff2.^(Nv(2)*(1-p2)).*nck1.*nck2));

% Entropy is given as the log of the number of states
entropy = log(omega);