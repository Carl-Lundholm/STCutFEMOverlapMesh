function value = PHI_km1k(a, b, K)

% Computes the local mass matrix element (phi_k, phi_km1) for a potentially 
% cut simplex K, i.e., with arbitrary integration interval (a, b), where
% x_km1 <= a <= b <= x_k.

x_km1 = K(1);
x_k = K(2);
h = K(3);

value = 1/(h^2)*((a^3 - b^3)/3 + (x_k + x_km1)*(b^2 - a^2)/2 + ...
    x_km1*x_k*(a - b));

