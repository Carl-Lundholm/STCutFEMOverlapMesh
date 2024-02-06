function b_fint = fIntegralTermdG1(K0, KG, x0, G, ...
                            amin, amax, bmin, bmax, mutr, t1ga, t1gb)

% Quadrature is used to approximate the integrals locally over 
% space-time prisms.
% 1st, temporal quadrature: composite 3-point Lobatto rule => quad error ~ k^4.
% 2nd, spatial quadrature: trapezoidal rule => quad error ~ h^2.

global kn x0_init x0_fin I0 tnm1 tn M leA

b_fint = zeros(leA, 1);

% Regular domains of U1 ---------------------------------------------------

tm = 0.5*(tnm1 + tn);
t_vec = [tnm1, tm, tn];

Crs = [1 2 0; 0 2 1];

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
    
    l = 1;
    bk1 = zeros(2,1);
    bk2 = zeros(2,1);
    for tl = t_vec
        
        F = f_func(alpha,tl) + f_func(beta,tl); 
        Fx = f_func(alpha,tl)*alpha + f_func(beta,tl)*beta;
        
        bk = [x_k*F - Fx; Fx - x_km1*F];    
    
        bk1 = bk1 + Crs(1,l)*bk;
        bk2 = bk2 + Crs(2,l)*bk;
        
        l = l + 1;
        
    end
    
    fac = (beta - alpha)*kn/(12*hk);
    
    bk1 = fac*bk1;
    bk2 = fac*bk2;
    
    b_fint(kpos) = b_fint(kpos) + bk1;
    b_fint(M + kpos) = b_fint(M + kpos) + bk2;
    
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
    
    alpha = max(x_km1, bmax);
    beta = x_k;
    
    l = 1;
    bk1 = zeros(2,1);
    bk2 = zeros(2,1);
    for tl = t_vec
        
        F = f_func(alpha,tl) + f_func(beta,tl); 
        Fx = f_func(alpha,tl)*alpha + f_func(beta,tl)*beta;
        
        bk = [x_k*F - Fx; Fx - x_km1*F];    
    
        bk1 = bk1 + Crs(1,l)*bk;
        bk2 = bk2 + Crs(2,l)*bk;
        
        l = l + 1;
        
    end
    
    fac = (beta - alpha)*kn/(12*hk);
    
    bk1 = fac*bk1;
    bk2 = fac*bk2;
    
    b_fint(kpos) = b_fint(kpos) + bk1;
    b_fint(M + kpos) = b_fint(M + kpos) + bk2;
    
end

% Gamma domains of U1 -----------------------------------------------------

% To the left of G
K1ga =  K0(:, K0(1,:) < amax & amin < K0(2,:));
leK1ga = length(K1ga(1,:));
for k = 1:leK1ga
    
    x_km1 = K1ga(1,k);
    x_k = K1ga(2,k);
    hk = K1ga(3, k);
    
    x_kpos = find((x0 == x_k));
    kpos = (x_kpos - 1: x_kpos);
    
    alpha = max(x_km1, amin);
    beta = min(x_k, amax);
    
    ti = t1ga(k);
    tf = t1ga(k+1);
    te = t1ga(end);
    tgm = 0.5*(ti + tf); 
    tm = 0.5*(tf + te);
    
    t_vec = [tf, tm, te];
    
    delt = abs(te - tf);
    deltg = abs(tf - ti);
    
    wgm = 2*deltg/3;
    wgf = deltg/6;
    wf = delt/6; 
    wm = 2*delt/3;
    we = wf;
    
    lamb1gm = tn - tgm;
    lamb1f = tn - tf;
    lamb1m = tn - tm;
    lamb1e = tn - te;
    
    lamb2gm = tgm - tnm1;
    lamb2f = tf - tnm1;
    lamb2m = tm - tnm1;
    lamb2e = te - tnm1;
    
    C1gm = wgm*lamb1gm;
    C1 = [(wgf + wf)*lamb1f, wm*lamb1m, we*lamb1e]; 
    
    C2gm = wgm*lamb2gm;
    C2 = [(wgf + wf)*lamb2f, wm*lamb2m, we*lamb2e]; 
        
    l = 1;
    bk1 = zeros(2,1);
    bk2 = zeros(2,1);
    for tl = t_vec
        
        F = f_func(alpha,tl) + f_func(beta,tl); 
        Fx = f_func(alpha,tl)*alpha + f_func(beta,tl)*beta;
        
        bk = [x_k*F - Fx; Fx - x_km1*F];    
    
        bk1 = bk1 + C1(l)*bk;
        bk2 = bk2 + C2(l)*bk;
        
        l = l + 1;
        
    end
    
    fac = (beta - alpha);
    
    bk1 = fac*bk1;
    bk2 = fac*bk2; 
    
    beta = 0.5*(alpha + beta);
    
    F = f_func(alpha,tgm) + f_func(beta,tgm); 
    Fx = f_func(alpha,tgm)*alpha + f_func(beta,tgm)*beta;
    
    bk = [x_k*F - Fx; Fx - x_km1*F];
    
    bk1 = bk1 + (beta - alpha)*C1gm*bk;
    bk2 = bk2 + (beta - alpha)*C2gm*bk;
    
    fac = 0.5/(kn*hk);  
    bk1 = fac*bk1;
    bk2 = fac*bk2;
    
    b_fint(kpos) = b_fint(kpos) + bk1;
    b_fint(M + kpos) = b_fint(M + kpos) + bk2;
    
    
end

% To the right of G
K1gb = K0(:, K0(1,:) < bmax & bmin < K0(2,:));
leK1gb = length(K1gb(1,:));
for k = 1:leK1gb
    
    x_km1 = K1gb(1,k);
    x_k = K1gb(2,k);
    hk = K1gb(3, k);
    
    x_kpos = find((x0 == x_k));
    kpos = (x_kpos - 1: x_kpos);
    
    alpha = max(x_km1, bmin);
    beta = min(x_k, bmax);
    
    ti = t1gb(k+1);
    tf = t1gb(k);
    te = t1gb(1);
    tgm = 0.5*(ti + tf); 
    tm = 0.5*(tf + te);
    
    t_vec = [tf, tm, te];
    
    delt = abs(te - tf);
    deltg = abs(tf - ti);
    
    wgm = 2*deltg/3;
    wgf = deltg/6;
    wf = delt/6; 
    wm = 2*delt/3;
    we = wf;
    
    lamb1gm = tn - tgm;
    lamb1f = tn - tf;
    lamb1m = tn - tm;
    lamb1e = tn - te;
    
    lamb2gm = tgm - tnm1;
    lamb2f = tf - tnm1;
    lamb2m = tm - tnm1;
    lamb2e = te - tnm1;
    
    C1gm = wgm*lamb1gm;
    C1 = [(wgf + wf)*lamb1f, wm*lamb1m, we*lamb1e]; 
    
    C2gm = wgm*lamb2gm;
    C2 = [(wgf + wf)*lamb2f, wm*lamb2m, we*lamb2e]; 
        
    l = 1;
    bk1 = zeros(2,1);
    bk2 = zeros(2,1);
    for tl = t_vec
        
        F = f_func(alpha,tl) + f_func(beta,tl); 
        Fx = f_func(alpha,tl)*alpha + f_func(beta,tl)*beta;
        
        bk = [x_k*F - Fx; Fx - x_km1*F];    
    
        bk1 = bk1 + C1(l)*bk;
        bk2 = bk2 + C2(l)*bk;
        
        l = l + 1;
        
    end
    
    fac = (beta - alpha);
    
    bk1 = fac*bk1;
    bk2 = fac*bk2; 
    
    alpha = 0.5*(alpha + beta);
    
    F = f_func(alpha,tgm) + f_func(beta,tgm); 
    Fx = f_func(alpha,tgm)*alpha + f_func(beta,tgm)*beta;
    
    bk = [x_k*F - Fx; Fx - x_km1*F];
    
    bk1 = bk1 + (beta - alpha)*C1gm*bk;
    bk2 = bk2 + (beta - alpha)*C2gm*bk;
    
    fac = 0.5/(kn*hk);  
    bk1 = fac*bk1;
    bk2 = fac*bk2;
    
    b_fint(kpos) = b_fint(kpos) + bk1;
    b_fint(M + kpos) = b_fint(M + kpos) + bk2;
end

% Domains of U2 -----------------------------------------------------------

tm = 0.5*(tnm1 + tn);
t_vec = [tnm1, tm, tn];

Crs = [1 2 0; 0 2 1];

K2r = KG(:, x0_init <= KG(1,:) & KG(2,:) <= x0_fin);
leK2r = length(K2r(1,:));
for k = 1:leK2r
    
    x_km1 = K2r(1,k);
    x_k = K2r(2,k);
    
    x_kpos = find((G == x_k));
    kpos = I0 + (x_kpos - 1: x_kpos);
    
    alpha = x_km1 - mutr*(tn - t_vec);
    beta = x_k - mutr*(tn - t_vec);
    
    l = 1;
    bk1 = zeros(2,1);
    bk2 = zeros(2,1);
    for tl = t_vec
        
        F = f_func(alpha(l),tl) + f_func(beta(l),tl); 
        Fx = f_func(alpha(l),tl)*alpha(l) + f_func(beta(l),tl)*beta(l);
        
        bk = [beta(l)*F - Fx; Fx - alpha(l)*F];    
    
        bk1 = bk1 + Crs(1,l)*bk;
        bk2 = bk2 + Crs(2,l)*bk;
        
        l = l + 1;
        
    end
    
    fac = kn/12;
    
    bk1 = fac*bk1;
    bk2 = fac*bk2;
    
   b_fint(kpos) = b_fint(kpos) + bk1;
   b_fint(M + kpos) = b_fint(M + kpos) + bk2;
    
end
    
    
    
    
    
    
    
    
    
    