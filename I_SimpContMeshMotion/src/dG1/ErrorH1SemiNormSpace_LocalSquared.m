function grad_error = ErrorH1SemiNormSpace_LocalSquared(grad_U, u, a, b)

% This function computes the local squared H1-seminorm of the error on the 
% intra-simplicial interval (a, b).
% The approximative function gradient grad_U should be given as a constant.
% The exact funtion u should be given as a Matlab-function that returns
% three values where the second value is the gradient of the function u.
% The integral is computed with 3-point Gauss-Legendre quadrature.

% Quadrature points for reference interval [-1, 1]
quad_points_ref = [-sqrt(3/5); 0; sqrt(3/5)];

% Quadrature weights
quad_weights = [5/9; 8/9; 5/9];

% Initialize error, actually error^2
grad_error = 0;

% Compute the integral of (u - U)^2 from a to b with 3-point 
% Gauss-Legendre quadrature 
delx = b - a;
x_mid = 0.5*(b + a);

for i = 1:3
    x_quad = x_mid + 0.5*delx*quad_points_ref(i);
    [~, grad_u_i, ~] = u(x_quad);    
    grad_error_i = (grad_u_i - grad_U)^2;
    grad_error = grad_error + grad_error_i*quad_weights(i);      
end
    
grad_error = 0.5*delx*grad_error;

    


