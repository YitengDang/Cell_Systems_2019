function [omega, omegak] = entropy_eq_sphere_2(dist_vec, Son, K, a0, Rcell)
% Calculate the analytical entropy given the parameters of the system
% dist_vec is a vector with all the distances r_i. dist_vec has N elements.
% omega k gives the number of eq. states for a given fraction of ON cells

dist_vec = dist_vec*a0;
N = numel(dist_vec);
r = dist_vec(dist_vec>0); % exclude self influence
fN = sum(sinh(Rcell)*sum(exp(Rcell-r)./r)); % calculate signaling strength
gN = sum(sum((sinh(Rcell)^2)*exp(2*(Rcell-r))./(r.^2))); % calculate noise variance strength

% nchoosek calculation is saved in data to reduce computation time. Then it
% is only calculated if the file does not exists yet. It saves the values
% of all Binom(N,n) for a given N and n = 0:N. The values are then used in
% the vector nck.
filename = fullfile(pwd,'data','nchoosek',sprintf('%d.mat', N));
if exist(filename, 'file') == 2
    tmp = load(filename);
    nck = tmp.nck;
else
    nck = zeros(N+1, 1);
    for n = 0:N
        nck(n+1) = nchoosek(N,n);
    end
    save(filename, 'nck');
end

% Calculate for all fractions
p = (0:N)/N;
muON = Son + fN*(p*Son+(1-p));
muOFF = 1 + fN*(p*Son+(1-p));
sigmap = sqrt(p.*(1-p)*gN)*(Son-1);
omegak = ((normcdf((K-muOFF)./sigmap)).^((1-p)*N) ...
    .* (normcdf((-K+muON)./sigmap)).^(p*N)).*nck';
omega = sum(sum(omegak));

% Entropy is given as the log of the number of states
omega = log(omega);