function error = ErrorXNormSlabwiseSquared_AuxiliaryOverlap(Un, Unm1p)

% Quadrature is used to approximate the integrals locally over 
% space-time prisms.
% 1st, temporal quadrature: 3-point Gauss-Legendre rule
% => quad error ~ (k^6)^(1/2) = k^3.
% 2nd, NO spatial quadrature! Exact computation since [grad e] = -[grad U]

global K0 KG amin amax bmin bmax t1ga t1gb I0 M tnm1 tn x0 mutr

% Quadrature points for reference interval [-1, 1]
quad_points_ref = [-sqrt(3/5); 0; sqrt(3/5)];

% Quadrature weights
quad_weights = [5/9; 8/9; 5/9];

% Initialize error
error = 0;

% Left overlap part
Kgs = K0(:, K0(1,:) < amax & amin < K0(2,:));
gmin = amin;
t1gg = t1ga;
indGL = I0 + 1;
indGR = I0 + 2;
hGb = KG(3, 1);

kg = 0; % Local index of Gamma-cut simplex in set Ks
for K = Kgs
    
    kg = kg + 1;
    
    x_km1 = K(1);
    x_k = K(2);
    hk = K(3);
    
    k = find((x0 == x_k));
    
    alpha0 = max(x_km1, gmin);
    beta0 = x_k;    
    
    tmid = 0.5*(t1gg(kg) + t1gg(kg+1));
    delt = abs(t1gg(kg+1) - t1gg(kg));
    
    quad_K = 0;
    for m = 1:3
       
        t_quad = tmid + 0.5*delt*quad_points_ref(m);
        
        % Compute (grad(U))_1
        U1_km1 = LinearInterpolationSingVal(t_quad, tnm1, tn, ...
                                        Unm1p(k-1), Un(k-1));
        U1_k = LinearInterpolationSingVal(t_quad, tnm1, tn, ...
                                        Unm1p(k), Un(k));
        grad_U1 = (U1_k - U1_km1)/hk;
        
        % Compute (grad(U))_2
        U2L = LinearInterpolationSingVal(t_quad, tnm1, tn, ...
                                        Unm1p(indGL), Un(indGL));
        U2R = LinearInterpolationSingVal(t_quad, tnm1, tn, ...
                                        Unm1p(indGR), Un(indGR));
        grad_U2 = (U2R - U2L)/hGb;
        
        % Compute length of overlap interval (SPECIAL FOR LEFT PART)
        alpha = alpha0 + mutr*(t_quad - t1gg(kg));
        length_Overlap_Km = beta0 - alpha;
   
        % Compute current quad term and add to quad for simplex K
        quad_Km = (grad_U1 - grad_U2)^2*length_Overlap_Km;
        quad_K = quad_K + quad_Km*quad_weights(m);
        
    end
    
    error = error + 0.5*delt*quad_K; 
        
end

% Right overlap part
Kgs = K0(:, K0(1,:) < bmax & bmin < K0(2,:));
gmin = bmin;
t1gg = t1gb;
indGL = M - 1;
indGR = M;
hGb = KG(3, end);

kg = 0; % Local index of Gamma-cut simplex in set Ks
for K = Kgs
    
    kg = kg + 1;
    
    x_km1 = K(1);
    x_k = K(2);
    hk = K(3);
    
    k = find((x0 == x_k));
    
    alpha0 = x_km1;
    beta0 = max(x_km1, gmin);    
    
    tmid = 0.5*(t1gg(kg) + t1gg(kg+1));
    delt = abs(t1gg(kg+1) - t1gg(kg));
    
    quad_K = 0;
    for m = 1:3
       
        t_quad = tmid + 0.5*delt*quad_points_ref(m);
        
        % Compute (grad(U))_1
        U1_km1 = LinearInterpolationSingVal(t_quad, tnm1, tn, ...
                                        Unm1p(k-1), Un(k-1));
        U1_k = LinearInterpolationSingVal(t_quad, tnm1, tn, ...
                                        Unm1p(k), Un(k));
        grad_U1 = (U1_k - U1_km1)/hk;
        
        % Compute (grad(U))_2
        U2L = LinearInterpolationSingVal(t_quad, tnm1, tn, ...
                                        Unm1p(indGL), Un(indGL));
        U2R = LinearInterpolationSingVal(t_quad, tnm1, tn, ...
                                        Unm1p(indGR), Un(indGR));
        grad_U2 = (U2R - U2L)/hGb;
        
        % Compute length of overlap interval (SPECIAL FOR RIGHT PART)
        beta = beta0 + mutr*(t_quad - t1gg(kg));
        length_Overlap_Km = beta - alpha0;
   
        % Compute current quad term and add to quad for simplex K
        quad_Km = (grad_U1 - grad_U2)^2*length_Overlap_Km;
        quad_K = quad_K + quad_Km*quad_weights(m);
        
    end
    
    error = error + 0.5*delt*quad_K; 
        
end



end

