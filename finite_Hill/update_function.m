function x_new = update_function(x, n, Con, K)
x_new = ((Con-1)*x + 1)^n/(K^n+((Con-1)*x + 1)^n) - x;