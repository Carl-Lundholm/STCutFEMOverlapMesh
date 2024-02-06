function value = PHI_km1lm1(a, b, K_k, K_l)

% Computes the local mass matrix element (phi_lm1, phi_km1) for 
% potentially cut simplices K_k and K_l, i.e., with arbitrary integration 
% interval (a, b), where x_lm1 <= a <= b <= x_l.

x_k = K_k(2);
h_k = K_k(3);

x_l = K_l(2);
h_l = K_l(3);

value = 1/(h_l*h_k)*((b^3 - a^3)/3 + 0.5*(x_l + x_k)*(a^2 - b^2) + ...
    x_l*x_k*(b - a));

