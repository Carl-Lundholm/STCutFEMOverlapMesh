function quad = Quad3pGL(f, a, b)

% This function approximates the integral of "f" from "a" to "b" by
% using 3-point Gauss-Legendre quadrature

% Quadrature points for reference interval [-1, 1]
quad_points_ref = [-sqrt(3/5); 0; sqrt(3/5)];

% Quadrature weights
quad_weights = [5/9; 8/9; 5/9];

% Initialize quadrature sum
quad = 0;

% Compute the integral of f from a to b with 3-point Gauss-Legendre quad 
delx = b - a;
x_mid = 0.5*(b + a);

for i = 1:3
    x_i = x_mid + 0.5*delx*quad_points_ref(i);
    quad = quad + f(x_i)*quad_weights(i);      
end
    
quad = 0.5*delx*quad;