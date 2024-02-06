function value = PHI_klm1(a, b, K_k, K_l)

% Computes the local mass matrix element (phi_lm1, phi_k) for 
% potentially cut simplices K_k and K_l, i.e., with arbitrary integration 
% interval (a, b), where x_lm1 <= a <= b <= x_l.

x_km1 = K_k(1);
h_k = K_k(3);

x_l = K_l(2);
h_l = K_l(3);

value = 1/(h_l*h_k)*((a^3 - b^3)/3 + 0.5*(x_km1 + x_l)*(b^2 - a^2) + ...
    x_km1*x_l*(a - b));

