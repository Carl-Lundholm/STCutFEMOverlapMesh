function error = ErrorXNormSlabwiseSquared_AuxiliaryInterior(Un, Unm1p, u)

% Quadrature is used to approximate the integrals locally over 
% space-time prisms.
% 1st, temporal quadrature: 3-point Gauss-Legendre rule
% => quad error ~ (k^6)^(1/2) = k^3.
% 2nd, spatial quadrature: 3-point Gauss-Legendre rule
% => quad error ~ (h^6)^(1/2) = h^3.

global x0_init x0_fin I0 tnm1 tn kn x0 K0 G KG amin amax bmin bmax mutr ...
    t1ga t1gb kn1ga kn1gb

% Quadrature points for reference interval [-1, 1]
quad_points_ref = [-sqrt(3/5); 0; sqrt(3/5)];

% Quadrature weights
quad_weights = [5/9; 8/9; 5/9];

error = 0;

% Uncut cells of Omega_1 --------------------------------------------------

tmid = 0.5*(tnm1 + tn);
delt = kn;

% To the left of G
Ks = K0(:, K0(1,:) < amin);
for K = Ks
    
    x_km1 = K(1);
    x_k = K(2);
    hk = K(3);
    
    k = find((x0 == x_k));
    
    alpha = x_km1;
    beta = min(x_k, amin);
    
    % D_t U
    DtU_km1 = (Un(k-1) - Unm1p(k-1))/kn;
    DtU_k = (Un(k) - Unm1p(k))/kn;
    DtU_K = [DtU_km1; DtU_k];
    
    quad_K = 0;
    for m = 1:3
       
        t_quad = tmid + 0.5*delt*quad_points_ref(m);
  
        % |D_t error|^2
        quad_Dt_Km = ErrorDtSpace_LocalSquared(DtU_K, K, @(x)u(x, t_quad), ... 
                                                alpha, beta, 0);
        
        % |grad error|^2
        U_km1 = LinearInterpolationSingVal(t_quad, tnm1, tn, ...
                                                  Unm1p(k-1), Un(k-1));
        U_k = LinearInterpolationSingVal(t_quad, tnm1, tn, Unm1p(k), Un(k));     
        grad_U = (U_k - U_km1)/hk;
        
        quad_grad_Km = ErrorH1SemiNormSpace_LocalSquared(grad_U, @(x)u(x, t_quad), ... 
                                                alpha, beta);
                                            
        % Add and update                                
        quad_K = quad_K + (kn*quad_Dt_Km + quad_grad_Km)*quad_weights(m);
        
    end
    
    error = error + 0.5*delt*quad_K;
    
end

% To the right of G
Ks = K0(:, bmax < K0(2,:));
for K = Ks
    
    x_km1 = K(1);
    x_k = K(2);
    hk = K(3);
    
    k = find((x0 == x_k));
    
    alpha = max(x_km1, bmax);
    beta = x_k;
    
    % D_t U
    DtU_km1 = (Un(k-1) - Unm1p(k-1))/kn;
    DtU_k = (Un(k) - Unm1p(k))/kn;
    DtU_K = [DtU_km1; DtU_k];
    
    quad_K = 0;
    for m = 1:3
       
        t_quad = tmid + 0.5*delt*quad_points_ref(m);
        
        % |D_t error|^2
        quad_Dt_Km = ErrorDtSpace_LocalSquared(DtU_K, K, @(x)u(x, t_quad), ... 
                                                alpha, beta, 0);
        
        % |grad error|^2
        U_km1 = LinearInterpolationSingVal(t_quad, tnm1, tn, ...
                                                  Unm1p(k-1), Un(k-1));
        U_k = LinearInterpolationSingVal(t_quad, tnm1, tn, Unm1p(k), Un(k));     
        grad_U = (U_k - U_km1)/hk;
        
        quad_grad_Km = ErrorH1SemiNormSpace_LocalSquared(grad_U, @(x)u(x, t_quad), ... 
                                                alpha, beta);
                                            
        % Add and update                                
        quad_K = quad_K + (kn*quad_Dt_Km + quad_grad_Km)*quad_weights(m);                                  
        
    end
    
    error = error + 0.5*delt*quad_K;
    
end

% Cut cells of Omega_1 ----------------------------------------------------

% To the left of G
Ks = K0(:, K0(1,:) < amax & amin < K0(2,:));
kg = 0; % Local index of Gamma-cut simplex in set Ks
for K = Ks
    
    kg = kg + 1;
    
    x_km1 = K(1);
    x_k = K(2);
    hk = K(3);
    
    k = find((x0 == x_k));
    
    alpha0 = max(x_km1, amin);
    beta0 = min(x_k, amax);
    
    % D_t U
    DtU_km1 = (Un(k-1) - Unm1p(k-1))/kn;
    DtU_k = (Un(k) - Unm1p(k))/kn;
    DtU_K = [DtU_km1; DtU_k];
    
    % Triangular part of space-time prism
    tmid = 0.5*(t1ga(kg) + t1ga(kg+1));
    delt = abs(t1ga(kg+1) - t1ga(kg));
    quad_K = 0;
    for m = 1:3
       
        t_quad = tmid + 0.5*delt*quad_points_ref(m);
        beta = alpha0 + mutr*(t_quad - t1ga(kg));
        
        % |D_t error|^2
        quad_Dt_Km = ErrorDtSpace_LocalSquared(DtU_K, K, @(x)u(x, t_quad), ... 
                                                alpha0, beta, 0);
        
        % |grad error|^2
        U_km1 = LinearInterpolationSingVal(t_quad, tnm1, tn, ...
                                                  Unm1p(k-1), Un(k-1));
        U_k = LinearInterpolationSingVal(t_quad, tnm1, tn, Unm1p(k), Un(k));     
        grad_U = (U_k - U_km1)/hk;
        
        quad_grad_Km = ErrorH1SemiNormSpace_LocalSquared(grad_U, @(x)u(x, t_quad), ... 
                                                alpha0, beta);
        
        % Add and update                                
        quad_K = quad_K + (kn*quad_Dt_Km + quad_grad_Km)*quad_weights(m);
        
    end
    
    error = error + 0.5*delt*quad_K;
    
    % Rectangular part of space-time prism
    tmid = kn1ga(2, kg + 1);
    delt = kn1ga(1, kg + 1);
    quad_K = 0;
    for m = 1:3
       
        t_quad = tmid + 0.5*delt*quad_points_ref(m);
        
        % |D_t error|^2
        quad_Dt_Km = ErrorDtSpace_LocalSquared(DtU_K, K, @(x)u(x, t_quad), ... 
                                                alpha0, beta0, 0);
        
        % |grad error|^2
        U_km1 = LinearInterpolationSingVal(t_quad, tnm1, tn, ...
                                                  Unm1p(k-1), Un(k-1));
        U_k = LinearInterpolationSingVal(t_quad, tnm1, tn, Unm1p(k), Un(k));     
        grad_U = (U_k - U_km1)/hk;
        
        quad_grad_Km = ErrorH1SemiNormSpace_LocalSquared(grad_U, @(x)u(x, t_quad), ... 
                                                alpha0, beta0);
                                            
        % Add and update                                
        quad_K = quad_K + (kn*quad_Dt_Km + quad_grad_Km)*quad_weights(m);
        
    end
    
    error = error + 0.5*delt*quad_K;
      
end

% To the right of G
Ks = K0(:, K0(1,:) < bmax & bmin < K0(2,:));
kg = 0; % Local index of Gamma-cut simplex in set Ks
for K = Ks
    
    kg = kg + 1;
    
    x_km1 = K(1);
    x_k = K(2);
    hk = K(3);
    
    k = find((x0 == x_k));
    
    alpha0 = max(x_km1, bmin);
    beta0 = min(x_k, bmax);
    
    % D_t U
    DtU_km1 = (Un(k-1) - Unm1p(k-1))/kn;
    DtU_k = (Un(k) - Unm1p(k))/kn;
    DtU_K = [DtU_km1; DtU_k];
    
    % Triangular part of space-time prism
    tmid = 0.5*(t1gb(kg) + t1gb(kg+1));
    delt = abs(t1gb(kg + 1) - t1gb(kg));
    quad_K = 0;
    for m = 1:3
       
        t_quad = tmid + 0.5*delt*quad_points_ref(m);
        alpha = alpha0 + mutr*(t_quad - t1gb(kg));
  
        % |D_t error|^2
        quad_Dt_Km = ErrorDtSpace_LocalSquared(DtU_K, K, @(x)u(x, t_quad), ... 
                                                alpha, beta0, 0);
        
        % |grad error|^2      
        U_km1 = LinearInterpolationSingVal(t_quad, tnm1, tn, ...
                                                  Unm1p(k-1), Un(k-1));
        U_k = LinearInterpolationSingVal(t_quad, tnm1, tn, Unm1p(k), Un(k));     
        grad_U = (U_k - U_km1)/hk;
        
        quad_grad_Km = ErrorH1SemiNormSpace_LocalSquared(grad_U, @(x)u(x, t_quad), ... 
                                                alpha, beta0);
        
        % Add and update                                
        quad_K = quad_K + (kn*quad_Dt_Km + quad_grad_Km)*quad_weights(m);
        
    end
    
    error = error + 0.5*delt*quad_K;
    
    % Rectangular part of space-time prism
    tmid = kn1gb(2, kg);
    delt = kn1gb(1, kg);
    quad_K = 0;
    for m = 1:3
       
        t_quad = tmid + 0.5*delt*quad_points_ref(m);
        
        % |D_t error|^2
        quad_Dt_Km = ErrorDtSpace_LocalSquared(DtU_K, K, @(x)u(x, t_quad), ... 
                                                alpha0, beta0, 0);
        
        % |grad error|^2   
        U_km1 = LinearInterpolationSingVal(t_quad, tnm1, tn, ...
                                                  Unm1p(k-1), Un(k-1));
        U_k = LinearInterpolationSingVal(t_quad, tnm1, tn, Unm1p(k), Un(k));     
        grad_U = (U_k - U_km1)/hk;
        
        quad_grad_Km = ErrorH1SemiNormSpace_LocalSquared(grad_U, @(x)u(x, t_quad), ... 
                                                alpha0, beta0);
        
        % Add and update                                
        quad_K = quad_K + (kn*quad_Dt_Km + quad_grad_Km)*quad_weights(m);
        
    end
    
    error = error + 0.5*delt*quad_K;
    
end

% Cells of Omega_2 --------------------------------------------------------

tmid = 0.5*(tnm1 + tn);
delt = kn;

Ks = KG(:, x0_init <= KG(1,:) & KG(2,:) <= x0_fin);
for K = Ks
    
    x_km1 = K(1);
    x_k = K(2);
    hk = K(3);
    
    k = find((G == x_k));
    k = k + I0;
    
    % D_t U
    DtU_km1 = (Un(k-1) - Unm1p(k-1))/kn;
    DtU_k = (Un(k) - Unm1p(k))/kn;
    DtU_K = [DtU_km1; DtU_k];
    
    quad_K = 0;
    for m = 1:3
       
        t_quad = tmid + 0.5*delt*quad_points_ref(m);
        alpha = x_km1 - mutr*(tn - t_quad);
        beta = x_k - mutr*(tn - t_quad);
        
        % |D_t error|^2
        quad_Dt_Km = ErrorDtSpace_LocalSquared(DtU_K, [alpha; beta], @(x)u(x, t_quad), ... 
                                                alpha, beta, 1);
        
        % |grad error|^2  
        U_km1 = LinearInterpolationSingVal(t_quad, tnm1, tn, ...
                                                  Unm1p(k-1), Un(k-1));
        U_k = LinearInterpolationSingVal(t_quad, tnm1, tn, Unm1p(k), Un(k));     
        grad_U = (U_k - U_km1)/hk;
        
        quad_grad_Km = ErrorH1SemiNormSpace_LocalSquared(grad_U, @(x)u(x, t_quad), ... 
                                                alpha, beta);
                                            
        % Add and update                                
        quad_K = quad_K + (kn*quad_Dt_Km + quad_grad_Km)*quad_weights(m);
        
    end
    
    error = error + 0.5*delt*quad_K;
    

end

