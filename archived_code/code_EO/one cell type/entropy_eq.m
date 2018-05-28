function omega = entropy_eq(K,Son,N,a0)
% Calculate the analytical entropy given the parameters of the system with
% point-like cells
% dist_vec is a vector with all the distances r_i. dist_vec has N elements.
% omega k gives the number of eq. states for a given fraction of ON cells

aux = dist_hex(ceil(sqrt(N)),a0);
r = aux(aux>0);
fN = sum(sum(exp(-r)));
gN = sum(sum(exp(-2*r)));
omega = 0;

for k = 0:N
    p = k/N;
    muON = Son + fN*(p*Son+(1-p));
    muOFF = 1 + fN*(p*Son+(1-p));
    sigmap = sqrt(p*(1-p)*gN*(Son-1)^2);
    omega = omega + (((1-p)*normcdf((K-muOFF)/sigmap) + ...
        p*(1-normcdf((K-muON)/sigmap)))^N)*nchoosek(N,k);
end

omega = omega/2^N;