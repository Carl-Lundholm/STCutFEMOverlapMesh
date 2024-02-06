function value = PHI_kk(a, b, K)

% Computes the local mass matrix element (phi_k, phi_k) for a potentially 
% cut simplex K, i.e., with arbitrary integration interval (a, b), where
% x_km1 <= a <= b <= x_k.

x_km1 = K(1);
h = K(3);

value = 1/(3*h^2)*((b - x_km1)^3 - (a - x_km1)^3);

