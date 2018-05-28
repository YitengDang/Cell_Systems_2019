function I = calc_I(p, theta, fN)
% Calculate the spatial order parameter given p and theta
eps = 1e-6;
if p < eps || p > 1-eps
    I = 0;
else
    I = (theta - (2*p-1)^2*fN)/4/p/(1-p)/fN;
end

    