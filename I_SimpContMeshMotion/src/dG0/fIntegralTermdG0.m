function b_fint = fIntegralTermdG0(K0, KG, x0, G, ...
                                amin, amax, bmin, bmax, mutr, t1ga, t1gb)

% Quadrature is used to approximate the integrals locally over 
% space-time prisms.
% 1st, temporal quadrature: midpoint rule => quad error ~ k^2.
% 2nd, spatial quadrature: trapezoidal rule => quad error ~ h^2.
                            
global kn kn1ga kn1gb x0_init x0_fin I0 tnm1 tn IG

b_fint = zeros(I0 + IG, 1);

t_av = 0.5*(tn + tnm1);         % Time that is used in the quadrature

% Regular domains of U1 ---------------------------------------------------

% To the left of G
K1ra = K0(:, K0(1,:) < amin);
leK1ra = length(K1ra(1,:));
for k = 1:leK1ra
    
    x_km1 = K1ra(1,k);
    x_k = K1ra(2,k);
    hk = K1ra(3, k);
    
    x_kpos = find((x0 == x_k));
    kpos = (x_kpos - 1: x_kpos);
    
    alpha = x_km1;
    beta = min(x_k, amin);
    
    delx = beta - alpha;
    
    alpha_kn = 0.5*(f_func(alpha,t_av) + f_func(beta,t_av));
    beta_kn = 0.5*(f_func(alpha,t_av)*alpha + f_func(beta,t_av)*beta);
    
    bk = 1/hk*delx*kn*[x_k*alpha_kn - beta_kn; beta_kn - x_km1*alpha_kn];
    
    b_fint(kpos) = b_fint(kpos) + bk;
    
end

% To the right of G
K1rb = K0(:, bmax < K0(2,:));
leK1rb = length(K1rb(1,:));
for k = 1:leK1rb
    
    x_km1 = K1rb(1,k);
    x_k = K1rb(2,k);
    hk = K1rb(3, k);
    
    x_kpos = find((x0 == x_k));
    kpos = (x_kpos - 1: x_kpos);
    
    alpha = max(bmax, x_km1);
    beta = x_k;
    
    delx = beta - alpha;
    
    alpha_kn = 0.5*(f_func(alpha,t_av) + f_func(beta,t_av));
    beta_kn = 0.5*(f_func(alpha,t_av)*alpha + f_func(beta,t_av)*beta);
      
    bk = 1/hk*delx*kn*[x_k*alpha_kn - beta_kn; beta_kn - x_km1*alpha_kn];
    
    b_fint(kpos) = b_fint(kpos) + bk;
    
end

% Gamma domains of U1 -----------------------------------------------------

% To the left of G
K1ga =  K0(:, amin < K0(2,:) & K0(1,:) < amax);         
leK1ga = length(K1ga(1,:));
for k = 1:leK1ga
    
    x_km1 = K1ga(1,k);
    x_k = K1ga(2,k);
    hk = K1ga(3, k);
    
    x_kpos = find((x0 == x_k));
    kpos = (x_kpos - 1: x_kpos);
    
    % Midpoint rule x 2. Both rectangular and triangular part
    ti = t1ga(k);
    tf = t1ga(k+1);
    te = t1ga(end);
    tgm = 0.5*(ti + tf); 
    tm = 0.5*(tf + te);
    
    delt = abs(te - tf);
    deltg = abs(tf - ti);
       
    % Rectangular part of space-time prism
    alpha = max(amin, x_km1);
    beta = min(amax, x_k);  
    delx = beta - alpha;
    
    F = f_func(alpha,tm) + f_func(beta,tm); 
    Fx = f_func(alpha,tm)*alpha + f_func(beta,tm)*beta;
    
    bk = 0.5/hk*delx*delt*[x_k*F - Fx; Fx - x_km1*F];
    
    % Add triangular part of space-time prism
    beta = 0.5*(alpha + beta);
    delx = beta - alpha;
    
    F = f_func(alpha,tgm) + f_func(beta,tgm); 
    Fx = f_func(alpha,tgm)*alpha + f_func(beta,tgm)*beta;
    
    bk = bk + 0.5/hk*delx*deltg*[x_k*F - Fx; Fx - x_km1*F];
     
%     % Weird mixed space-time "midpoint"-trapezoidal quadrature. 
%     % Not sure about quad error order. Worst case is k^1 + h^1.
%     intdom = 0.5*delx*(kn1ga(1,k) + kn1ga(1,k+1));
%  
%     t_alpha = kn1ga(2,k);
%     t_beta = kn1ga(2,k+1);
%     
%     alpha_kn = 0.5*(f_func(alpha,t_alpha) + f_func(beta,t_beta));
%     beta_kn = 0.5*(f_func(alpha,t_alpha)*alpha + f_func(beta,t_beta)*beta);
% 
%     bk = intdom/hk*[x_k*alpha_kn - beta_kn; beta_kn - x_km1*alpha_kn];
    
    % Add local contribution to global vector
    b_fint(kpos) = b_fint(kpos) + bk;
    
end

% To the right of G
K1gb =  K0(:, bmin < K0(2,:) & K0(1,:) < bmax);           
leK1gb = length(K1gb(1,:));
for k = 1:leK1gb
    
    x_km1 = K1gb(1,k);
    x_k = K1gb(2,k);
    hk = K1gb(3, k);
    
    x_kpos = find((x0 == x_k));
    kpos = (x_kpos - 1: x_kpos);
    
    % Midpoint rule x 2. Both rectangular and triangular part
    ti = t1gb(k+1);
    tf = t1gb(k);
    te = t1gb(1);
    tgm = 0.5*(ti + tf); 
    tm = 0.5*(tf + te);
    
    delt = abs(te - tf);
    deltg = abs(tf - ti);
       
    % Rectangular part of space-time prism
    alpha = max(bmin, x_km1);
    beta = min(bmax, x_k);
    delx = beta - alpha;
    
    F = f_func(alpha,tm) + f_func(beta,tm); 
    Fx = f_func(alpha,tm)*alpha + f_func(beta,tm)*beta;
    
    bk = 0.5/hk*delx*delt*[x_k*F - Fx; Fx - x_km1*F];
    
    % Add triangular part of space-time prism
    alpha = 0.5*(alpha + beta);
    delx = beta - alpha;
    
    F = f_func(alpha,tgm) + f_func(beta,tgm); 
    Fx = f_func(alpha,tgm)*alpha + f_func(beta,tgm)*beta;
    
    bk = bk + 0.5/hk*delx*deltg*[x_k*F - Fx; Fx - x_km1*F];
    
%     % Weird mixed space-time "midpoint"-trapezoidal quadrature. 
%     % Not sure about quad error order. Worst case is k^1 + h^1.
%     intdom = 0.5*delx*(kn1gb(1,k) + kn1gb(1,k+1));
%    
%     t_alpha = kn1gb(2,k);
%     t_beta = kn1gb(2,k+1);
%       
%     alpha_kn = 0.5*(f_func(alpha,t_alpha) + f_func(beta,t_beta));
%     beta_kn = 0.5*(f_func(alpha,t_alpha)*alpha + f_func(beta,t_beta)*beta);
%     
%     bk = intdom/hk*[x_k*alpha_kn - beta_kn; beta_kn - x_km1*alpha_kn];
    
    % Add local contribution to global vector
    b_fint(kpos) = b_fint(kpos) + bk;
    
end

% Domains of U2 -----------------------------------------------------------

K2r = KG(:, x0_init <= KG(1,:) & KG(2,:) <= x0_fin);
leK2r = length(K2r(1,:));
for k = 1:leK2r
    
    x_km1 = K2r(1,k);
    x_k = K2r(2,k);
    
    x_kpos = find((G == x_k));
    kpos = I0 + (x_kpos - 1: x_kpos);
    
    x_km1q = x_km1 - mutr*(tn - t_av);
    x_kq = x_k - mutr*(tn - t_av);
    
    alpha_kn = 0.5*(f_func(x_km1q,t_av) + f_func(x_kq,t_av));
    beta_kn = 0.5*(f_func(x_km1q,t_av)*x_km1q + f_func(x_kq,t_av)*x_kq);
    
    bk = kn*[x_kq*alpha_kn - beta_kn; beta_kn - x_km1q*alpha_kn];
    
    b_fint(kpos) = b_fint(kpos) + bk;
    
end


    
    
    
    
    
    
    
    
    
    