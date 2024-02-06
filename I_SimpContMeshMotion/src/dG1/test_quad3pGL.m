% Test the 3-point Gauss-Legendre quadrature
% Should be exact for 2(3) - 1 = 5 degree polynomials

% Define integration domain
a = 0; b = 1;

% Compute quadrature errors for random degree n polynomials
for n = 0:7
    
    % Construct a random degree n polynomial
    as = rand(n+1, 1);
    polynom = @(x)(0);
    for k = 0:n
        polynom = @(x)(polynom(x) + as(k+1)*x^(k));
    end
    
    % Compute exact integral
    integral_exact = 0;
    for k = 0:n
        term_n = as(k+1)/(k+1)*(b^(k+1) - a^(k+1));
        
        integral_exact = integral_exact + term_n;
    end
    
    % Compute approximate integral
    % NOTE: The input arguments of U and u have be modified since
    % the function below computes the squared L2-error, i.e., the integral
    % from a to b of (U - u)^2. With U = 0 and u = sqrt(polynom), the
    % integrand is "polynom" instead.
    integral_quad = ErrorL2NormSpaceQuad3PointGaussLegendre_Auxiliary( ...
        [a, b], [0, 0], @(x)(sqrt(polynom(x))));
    
    error_quad = integral_quad - integral_exact;
    
    disp(['Polynomial degree = ', num2str(n), ...
            ', Quad error = ', num2str(error_quad)])
      
end 



