function error = ErrorXNormSlabwiseSquared_AuxiliaryGammaGen(Un, Unm1p, ...
                 u, nt_abs, nx_abs, Kgs, gmin, gmax, t1gg, ...
                 indG, indGL, indGR, hGb)

% Quadrature is used to approximate the integrals locally over each
% Gamma-segment in a space-time prism.
% Spatio-temporal quadrature along a Gamma-segment: 
% 3-point Gauss-Legendre rule => quad error ~ ((dels)^6)^(1/2) = dels^3
% ((k^2 + h^2)^(1/2))^3 ~ k^3 + h^3

global tnm1 tn x0 C w1 w2

% Quadrature points for reference interval [-1, 1]
quad_points_ref = [-sqrt(3/5); 0; sqrt(3/5)];

% Quadrature weights
quad_weights = [5/9; 8/9; 5/9];

error = 0;

kg = 0; % Local index of Gamma-cut simplex in set Ks
for K = Kgs
    
    kg = kg + 1;
    
    x_km1 = K(1);
    x_k = K(2);
    hk = K(3);
    
    k = find((x0 == x_k));
    
    alpha = max(x_km1, gmin);
    beta = min(x_k, gmax);    
    xmid = 0.5*(alpha + beta);
    delx = beta - alpha;
    
    tmid = 0.5*(t1gg(kg) + t1gg(kg+1));
    delt = abs(t1gg(kg+1) - t1gg(kg));
    
    dels = C*delt;
    
    quad_K = 0;
    quadGrad_K = 0;
    for m = 1:3
       
        x_quad = xmid + 0.5*delx*quad_points_ref(m);
        t_quad = tmid + 0.5*delt*quad_points_ref(m);
         
        U1_km1 = LinearInterpolationSingVal(t_quad, tnm1, tn, ...
                                                  Unm1p(k-1), Un(k-1));
        U1_k = LinearInterpolationSingVal(t_quad, tnm1, tn, Unm1p(k), Un(k));
        
        % Quadrature for terms ||[e]||^2_{1/2,h,Gamma_n} and
        % ||abs(nt)^{1/2}[e]||^2_{Gamma_n} 
        U1 = LinearInterpolationSingVal(x_quad, x_km1, x_k, U1_km1, U1_k); 
        
        U2 = LinearInterpolationSingVal(t_quad, tnm1, tn, ...
                                        Unm1p(indG), Un(indG));
   
        quad_Km = (U1 - U2)^2;
        quad_K = quad_K + quad_Km*quad_weights(m);
        
        % Quadrature for term ||<dne>||^2_{-1/2,h,Gamma_n}        
        grad_U1 = (U1_k - U1_km1)/hk;
        U2L = LinearInterpolationSingVal(t_quad, tnm1, tn, ...
                                        Unm1p(indGL), Un(indGL));
        U2R = LinearInterpolationSingVal(t_quad, tnm1, tn, ...
                                        Unm1p(indGR), Un(indGR));
        grad_U2 = (U2R - U2L)/hGb;
        
        [~, grad_u] = u(x_quad, t_quad);
   
        quadGrad_Km = (grad_u - w1*grad_U1 - w2*grad_U2)^2;
        quadGrad_K = quadGrad_K + quadGrad_Km*quad_weights(m);
        
    end
    
    error = error + 0.5*dels*((1/hk + nt_abs)*quad_K + ...
                              nx_abs^2*hk*quadGrad_K); 
        
end

end

