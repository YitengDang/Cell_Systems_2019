function [omega, omegak] = entropy_eq_sphere(dist_vec, Son, K, a0, Rcell)
% Calculate the analytical entropy given the parameters of the system
% dist_vec is a vector with all the distances r_i. dist_vec has N elements.
% omega k gives the number of eq. states for a given fraction of ON cells

dist_vec = dist_vec*a0;
N = numel(dist_vec);
r = dist_vec(dist_vec>0); % exclude self influence
fN = sum(sinh(Rcell)*sum(exp(Rcell-r)./r)); % calculate signaling strength
gN = sum(sum((sinh(Rcell)^2)*exp(2*(Rcell-r))./(r.^2))); % calculate noise variance strength
omega = 0; % initialize the counting of number of states
omegak = zeros(N+1,1);
% Count the number of states given that k cells are ON
for k = 0:N
    p = k/N;
    muON = Son + fN*(p*Son+(1-p));
    muOFF = 1 + fN*(p*Son+(1-p));
    sigmap = sqrt(p*(1-p)*gN)*(Son-1);
%     omegak(k+1) = (((1-p)*normcdf((K-muOFF)/sigmap) + ...
%         p*(normcdf((-K+muON)/sigmap)))^N)*nchoosek(N,k);
    omegak(k+1) = ((normcdf((K-muOFF)/sigmap))^(N-k) ...
        * (normcdf((-K+muON)/sigmap))^k);
    omega = omega + omegak(k+1)*nchoosek(N,k);
end

% Normalize the entropy to the number of possible states
omega = log(omega);