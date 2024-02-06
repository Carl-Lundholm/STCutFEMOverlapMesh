function error = ErrorL2NormSpace(U, u, x0, G)

global I0 M

% This function computes the global spatial L2-norm of the error.
% The approximative function U should be given as a vector.
% The exact funtion u should be given as a MATLAB-function
% The integral is computed with composite 3-point Gauss-Legendre quadrature

% Endpoints of G
xGl = G(1);
xGr = G(end);

% Partition the function
U0 = U(1:I0)';
U2 = U(I0 + 1:M)';

% Initialize error, actually error^2
error = 0;

% Domains of U1 ---------------------------------------------------

% To the left of G
Xind = (x0 < xGl);
Xind(1) = 1;
xind_aR = find(Xind, 1, 'last' ) + 1;

xind_extr = xind_aR;
if xind_extr <= I0
    Xind(xind_extr) = 1;
end

xs = x0(Xind);
Us = U0(Xind);

Us(end) = LinearInterpolationSingVal(xGl, xs(end - 1), xs(end), ...
                                            Us(end - 1), Us(end));
xs(end) = xGl;

error_part = ErrorL2NormSpace_Auxiliary(Us, u, xs);

error = error + error_part;
    
% To the right of G
Xind = (x0 > xGr);
Xind(I0) = 1;
xind_bL = find(Xind, 1, 'first' ) - 1;

xind_extr = xind_bL;
if xind_extr > 0
    Xind(xind_extr) = 1;
end

xs = x0(Xind);
Us = U0(Xind);

Us(1) = LinearInterpolationSingVal(xGr, xs(1), xs(2), Us(1), Us(2));
xs(1) = xGr;

error_part = ErrorL2NormSpace_Auxiliary(Us, u, xs);

error = error + error_part;

% Domains of U2 -----------------------------------------------------------  

xs = G;
Us = U2;

error_part = ErrorL2NormSpace_Auxiliary(Us, u, xs);

error = error + error_part;

% Take the square root of the error 

error = sqrt(error);

% plot(G, err2, 'k-x');



