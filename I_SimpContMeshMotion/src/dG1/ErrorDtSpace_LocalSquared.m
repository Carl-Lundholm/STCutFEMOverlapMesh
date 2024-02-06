function Dt_error = ErrorDtSpace_LocalSquared(DtU_K, K, u, a, b, onG)

% This function computes the local squared L2-norm of the D_t of the error 
% on the intra-simplicial interval (a, b).
% The approximative function DtU_K should be given as a 2-vector containing
% the left and right endpoint values of DtU on K.
% The exact funtion u should be given as a Matlab-function that returns 
% three values where the second value is the gradient of u and the third 
% value is the time derivative of u
% The integral is computed with 3-point Gauss-Legendre quadrature.

global mutr

% Quadrature points for reference interval [-1, 1]
quad_points_ref = [-sqrt(3/5); 0; sqrt(3/5)];

% Quadrature weights
quad_weights = [5/9; 8/9; 5/9];

% Initialize error, actually error^2
Dt_error = 0;

% Compute the integral of (Dt(u - U))^2 from a to b with 3-point 
% Gauss-Legendre quadrature 
delx = b - a;
x_mid = 0.5*(b + a);

for i = 1:3
    x_quad = x_mid + 0.5*delx*quad_points_ref(i);
    [~, u_x, u_t] = u(x_quad);
    Dtu_i = u_t + onG*mutr*u_x;
    DtU_i = LinearInterpolationSingVal(x_quad, K(1), K(2), DtU_K(1), DtU_K(2));
    Dt_error_i = (Dtu_i - DtU_i)^2;
    Dt_error = Dt_error + Dt_error_i*quad_weights(i);      
end
    
Dt_error = 0.5*delx*Dt_error;

    


