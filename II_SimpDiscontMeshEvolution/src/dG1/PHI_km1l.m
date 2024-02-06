function value = PHI_km1l(a, b, K_k, K_l)

% Computes the local mass matrix element (phi_l, phi_km1) for 
% potentially cut simplices K_k and K_l, i.e., with arbitrary integration 
% interval (a, b), where x_lm1 <= a <= b <= x_l.

x_k = K_k(2);
h_k = K_k(3);

x_lm1 = K_l(1);
h_l = K_l(3);

value = 1/(h_l*h_k)*((a^3 - b^3)/3 + 0.5*(x_lm1 + x_k)*(b^2 - a^2) + ...
    x_lm1*x_k*(a - b));

