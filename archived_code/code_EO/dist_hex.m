function C0_X = dist_hex(n, a0)
% Old function from Theo's code to calculate the distances in hexagonal
% lattice

C0_X = zeros(n,n);
c = floor(n/2);
xc=sqrt(3)/2*c;
yc=(c/2 + c);
for k = 1 : n
    for j = 1 : n
        x = sqrt(3)/2*k;% coordinates transformation
        y = (k/2 + j);
        r = sqrt((x-xc)^2 + (y-yc)^2);
        C0_X(k,j) = r*a0;
    end
end