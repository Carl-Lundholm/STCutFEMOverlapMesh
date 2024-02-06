function value = PHI_km1km1(a, b, K)

% Computes the local mass matrix element (phi_km1, phi_km1) for a 
% potentially cut simplex K, i.e., with arbitrary integration interval 
% (a, b), where x_km1 <= a <= b <= x_k.

x_k = K(2);
h = K(3);

value = 1/(3*h^2)*((x_k - a)^3 - (x_k - b)^3);

