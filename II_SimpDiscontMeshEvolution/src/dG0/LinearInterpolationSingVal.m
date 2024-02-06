function y = LinearInterpolationSingVal(x, xp, xn, yp, yn)

y = (yn - yp)*(x - xp) *1/(xn - xp) + yp;  
