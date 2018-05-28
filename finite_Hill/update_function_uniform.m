function x_new = update_function_uniform(x, n, Con, K, fN)
x_new = (((1+fN)*(Con-1)*x + 1))^n/(K^n+(((1+fN)*(Con-1)*x + 1))^n) - x;
% for single cell in uniform lattice