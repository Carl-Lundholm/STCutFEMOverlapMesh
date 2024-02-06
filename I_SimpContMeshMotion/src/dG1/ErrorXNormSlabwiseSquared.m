function error = ErrorXNormSlabwiseSquared(Un, Unm1p, Unm1, u, n)

% This function computes the slabwise squared part of the X-norm of the error.
% The approximative function U should be given as a vector.
% The exact funtion u should be given as a MATLAB-function
% The integrals are computed with composite 3-point Gauss-Legendre quadrature
% in both space and time, which after taking the square root results in
% quad error ~ (k^6 + h^6)^(1/2) ~ k^3 + h^3.

global N tnm1 tn x0 G Gnm1

% Initialize the error of the current slab
error = 0;

% Compute slab terms of X-norm of the error

%% Terms involving integration over space time domains of codim = 0

term_eXn = ErrorXNormSlabwiseSquared_AuxiliaryInterior(Un, Unm1p, u);
error = error + term_eXn;

% We use regularity of u to write [grad(e)] = -[grad(U)] and thus 
% u is not supplied to the computation of the overlap terms.
term_eXn = ErrorXNormSlabwiseSquared_AuxiliaryOverlap(Un, Unm1p);
error = error + term_eXn;

%% Terms involving integration over space time domains of codim = 1

term_eXn = ErrorXNormSlabwiseSquared_AuxiliaryGamma(Un, Unm1p, u);
error = error + term_eXn;

%% Terms involving integration over pure spatial domains

if n == 2 

    % X-norm of error: Term ||e_0^+||^2
    term_eXn = ErrorL2NormSpace(Unm1p, @(x)(u(x,tnm1)), x0, Gnm1);  
    error = error + term_eXn^2;

else

    % X-norm of error: Term ||[e]_n-1||^2 
    % We use regularity of u to write [e] = -[U]
    % Thus letting the first argument be Unm1p - Unm1, and the second 0
    % in the function call below since this 
    term_eXn = ErrorL2NormSpace(Unm1p - Unm1, @(x)(0), x0, Gnm1);
    error = error + term_eXn^2;

end

if n == N
    
    % X-norm of error: Term ||e_N^-||^2
    term_eXn = ErrorL2NormSpace(Un, @(x)(u_func(x,tn)), x0, G);
    error = error + term_eXn^2;
    
end


