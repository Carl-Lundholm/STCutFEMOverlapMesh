function error = ErrorL2NormSpace_LocalSquared(Ua, Ub, u, a, b)

% This function computes the local squared L2-norm of the error on the 
% potentially intra-simplicial interval (a, b).
% The approximative function is defined by its endpoint values Ua and Ub
% The exact funtion u should be given as a Matlab-function
% The integral is computed with 3-point Gauss-Legendre quadrature

% Quadrature points for reference interval [-1, 1]
quad_points_ref = [-sqrt(3/5); 0; sqrt(3/5)];

% Quadrature weights
quad_weights = [5/9; 8/9; 5/9];

% Initialize error, actually error^2
error = 0;

% Compute the integral of (u - U)^2 from a to b with 3-point 
% Gauss-Legendre quadrature 
delx = b - a;
x_mid = 0.5*(b + a);

for i = 1:3
    x_quad = x_mid + 0.5*delx*quad_points_ref(i);
    U_ki = LinearInterpolationSingVal(x_quad, a, b, Ua, Ub);
    u_ki = u(x_quad);
    error_ki = abs(u_ki - U_ki)^2;
    error = error + error_ki*quad_weights(i);      
end
    
error = 0.5*delx*error;

    


