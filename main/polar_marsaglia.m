function x = polar_marsaglia
    s = 0;
    while s == 0 || s >=1
        tmp = 2*rand(2,1)-1;
        s = sum(tmp.^2);
    end
    x = tmp.*sqrt(-2*log(s)/s);