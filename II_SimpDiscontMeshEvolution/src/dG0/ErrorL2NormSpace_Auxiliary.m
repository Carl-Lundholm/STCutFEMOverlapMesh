function error = ErrorL2NormSpace_Auxiliary(U, u, x)

% This function computes the squared L2-norm of the error.
% The approximative function U should be given as a vector.
% The exact funtion u should be given as a Matlab-function
% The integral is computed with composite 3-point Gauss-Legendre quadrature

% Initialize error, actually error^2
error = 0;

% Compute the integral of (u - U)^2 w.r.t. x and over the xs with
% composite 3-point Gauss-Legendre quadrature 
for k = 2:length(x)
    
    local_eL2 = ErrorL2NormSpace_LocalSquared(U(k-1), U(k), u, x(k-1), x(k));
    
    error = error + local_eL2;
    
end
    


