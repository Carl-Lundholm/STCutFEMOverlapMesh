function b_fint = fIntegralTermdG0_cruderQuad(K0, KG, x0, G, ...
                                amin, amax, bmin, bmax)

global kn kn1ga kn1gb x0_init x0_fin I0 tnm1 tn IG mutr

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
    
    alpha_kn = 0.5*(f_func(x_km1,t_av) + f_func(x_k,t_av));
    beta_kn = 0.5*(f_func(x_km1,t_av)*x_km1 + f_func(x_k,t_av)*x_k);
    
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
    
    alpha_kn = 0.5*(f_func(x_km1,t_av) + f_func(x_k,t_av));
    beta_kn = 0.5*(f_func(x_km1,t_av)*x_km1 + f_func(x_k,t_av)*x_k);
      
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
    
    alpha = max(amin, x_km1);
    beta = min(amax, x_k);
    
    delx = beta - alpha;
    
    intdom = 0.5*delx*(kn1ga(1,k) + kn1ga(1,k+1));
 
    alpha_kn = 0.5*(f_func(x_km1,t_av) + f_func(x_k,t_av));
    beta_kn = 0.5*(f_func(x_km1,t_av)*x_km1 + f_func(x_k,t_av)*x_k);

    bk = intdom/hk*[x_k*alpha_kn - beta_kn; beta_kn - x_km1*alpha_kn];
    
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
    
    alpha = max(bmin, x_km1);
    beta = min(bmax, x_k);
    
    delx = beta - alpha;
    
    intdom = 0.5*delx*(kn1gb(1,k) + kn1gb(1,k+1));
   
    alpha_kn = 0.5*(f_func(x_km1,t_av) + f_func(x_k,t_av));
    beta_kn = 0.5*(f_func(x_km1,t_av)*x_km1 + f_func(x_k,t_av)*x_k);
    
    bk = intdom/hk*[x_k*alpha_kn - beta_kn; beta_kn - x_km1*alpha_kn];
    
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


    
    
    
    
    
    
    
    
    
    