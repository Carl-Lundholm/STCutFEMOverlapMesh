function A_dt = DtTermdG1(K0, KG, x0, G, ...
                                amin, amax, bmin, bmax, t1ga, t1gb)

global kn x0_init x0_fin I0 leA M tn tnm1

A_dt = zeros(leA);

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
    
    phi_km1km1 = -1/(3*hk^2)*((x_k - beta)^3 - (x_k - alpha)^3);
    phi_kkm1 = 1/hk^2*(-(beta^3 - alpha^3)/3 + ...
        (x_k + x_km1)*(beta^2 - alpha^2)/2 - x_km1*x_k*(beta - alpha));
    phi_kk = 1/(3*hk^2)*((beta - x_km1)^3 - (alpha - x_km1)^3);
    
    Ak = 0.5*[phi_km1km1, phi_kkm1; phi_kkm1, phi_kk];
      
    A_dt(kpos, kpos) = A_dt(kpos, kpos) - Ak;
    A_dt(kpos, M + kpos) = A_dt(kpos, M + kpos) + Ak;
    A_dt(M + kpos, kpos) = A_dt(M + kpos, kpos) - Ak;
    A_dt(M + kpos, M + kpos) = A_dt(M + kpos, M + kpos) + Ak;
   
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
    
    phi_km1km1 = -1/(3*hk^2)*((x_k - beta)^3 - (x_k - alpha)^3);
    phi_kkm1 = 1/hk^2*(-(beta^3 - alpha^3)/3 + ...
        (x_k + x_km1)*(beta^2 - alpha^2)/2 - x_km1*x_k*(beta - alpha));
    phi_kk = 1/(3*hk^2)*((beta - x_km1)^3 - (alpha - x_km1)^3);
    
    Ak = 0.5*[phi_km1km1, phi_kkm1; phi_kkm1, phi_kk];
      
    A_dt(kpos, kpos) = A_dt(kpos, kpos) - Ak;
    A_dt(kpos, M + kpos) = A_dt(kpos, M + kpos) + Ak;
    A_dt(M + kpos, kpos) = A_dt(M + kpos, kpos) - Ak;
    A_dt(M + kpos, M + kpos) = A_dt(M + kpos, M + kpos) + Ak;
    
end

% Gamma domains of U1 -----------------------------------------------------

% To the left of G
K1ga = K0(:, K0(1,:) < amax & amin < K0(2,:)); 
leK1ga = length(K1ga(1,:));
for k = 1:leK1ga
    
    x_km1 = K1ga(1,k);
    x_k = K1ga(2,k);
    hk = K1ga(3, k);
    
    x_kpos = find((x0 == x_k));
    kpos = (x_kpos - 1: x_kpos);
    
    alpha = max(x_km1, amin);
    beta = min(x_k, amax);
    
    phi_km1km1 = -1/(3*hk^2)*((x_k - beta)^3 - (x_k - alpha)^3);
    phi_kkm1 = 1/hk^2*(-(beta^3 - alpha^3)/3 + ...
        (x_k + x_km1)*(beta^2 - alpha^2)/2 - x_km1*x_k*(beta - alpha));
    phi_kk = 1/(3*hk^2)*((beta - x_km1)^3 - (alpha - x_km1)^3);
    
    phi_mat = [phi_km1km1, phi_kkm1; phi_kkm1, phi_kk];
    
    beta = 0.5*(alpha + beta);
    
    phi_km1km1gm = -1/(3*hk^2)*((x_k - beta)^3 - (x_k - alpha)^3);
    phi_kkm1gm = 1/hk^2*(-(beta^3 - alpha^3)/3 + ...
        (x_k + x_km1)*(beta^2 - alpha^2)/2 - x_km1*x_k*(beta - alpha));
    phi_kkgm = 1/(3*hk^2)*((beta - x_km1)^3 - (alpha - x_km1)^3);
    
    phi_matgm = [phi_km1km1gm, phi_kkm1gm; phi_kkm1gm, phi_kkgm];
    
    ti = t1ga(k);
    tf = t1ga(k+1);
    te = t1ga(end);
    tgm = 0.5*(ti + tf); 
    tm = 0.5*(tf + te);
    
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
    C1 = (wgf + wf)*lamb1f + wm*lamb1m + we*lamb1e; 
    
    C2gm = wgm*lamb2gm;
    C2 = (wgf + wf)*lamb2f + wm*lamb2m + we*lamb2e; 
    
    Ak1 = kn^-2*(C1gm*phi_matgm + C1*phi_mat);
    Ak2 = kn^-2*(C2gm*phi_matgm + C2*phi_mat);
      
    A_dt(kpos, kpos) = A_dt(kpos, kpos) - Ak1;
    A_dt(kpos, M + kpos) = A_dt(kpos, M + kpos) + Ak1;
    A_dt(M + kpos, kpos) = A_dt(M + kpos, kpos) - Ak2;
    A_dt(M + kpos, M + kpos) = A_dt(M + kpos, M + kpos) + Ak2;
    
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
    
    phi_km1km1 = -1/(3*hk^2)*((x_k - beta)^3 - (x_k - alpha)^3);
    phi_kkm1 = 1/hk^2*(-(beta^3 - alpha^3)/3 + ...
        (x_k + x_km1)*(beta^2 - alpha^2)/2 - x_km1*x_k*(beta - alpha));
    phi_kk = 1/(3*hk^2)*((beta - x_km1)^3 - (alpha - x_km1)^3);
    
    phi_mat = [phi_km1km1, phi_kkm1; phi_kkm1, phi_kk];
    
    alpha = 0.5*(alpha + beta);
    
    phi_km1km1gm = -1/(3*hk^2)*((x_k - beta)^3 - (x_k - alpha)^3);
    phi_kkm1gm = 1/hk^2*(-(beta^3 - alpha^3)/3 + ...
        (x_k + x_km1)*(beta^2 - alpha^2)/2 - x_km1*x_k*(beta - alpha));
    phi_kkgm = 1/(3*hk^2)*((beta - x_km1)^3 - (alpha - x_km1)^3);
    
    phi_matgm = [phi_km1km1gm, phi_kkm1gm; phi_kkm1gm, phi_kkgm];
    
    ti = t1gb(k+1);
    tf = t1gb(k);
    te = t1gb(1);
    tgm = 0.5*(ti + tf); 
    tm = 0.5*(tf + te);
    
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
    C1 = (wgf + wf)*lamb1f + wm*lamb1m + we*lamb1e; 
    
    C2gm = wgm*lamb2gm;
    C2 = (wgf + wf)*lamb2f + wm*lamb2m + we*lamb2e; 
    
    Ak1 = kn^-2*(C1gm*phi_matgm + C1*phi_mat);
    Ak2 = kn^-2*(C2gm*phi_matgm + C2*phi_mat);
      
    A_dt(kpos, kpos) = A_dt(kpos, kpos) - Ak1;
    A_dt(kpos, M + kpos) = A_dt(kpos, M + kpos) + Ak1;
    A_dt(M + kpos, kpos) = A_dt(M + kpos, kpos) - Ak2;
    A_dt(M + kpos, M + kpos) = A_dt(M + kpos, M + kpos) + Ak2;  

end

% Domain of U2 ------------------------------------------------------------

K2r = KG(:, x0_init <= KG(1,:) & KG(2,:) <= x0_fin);
leK2r = length(K2r(1,:));
for k = 1:leK2r
    
    x_km1 = K2r(1,k);
    x_k = K2r(2,k);
    hk = K2r(3, k);
    
    x_kpos = find((G == x_k));
    kpos = I0 + (x_kpos - 1: x_kpos);
    
    alpha = x_km1;
    beta = x_k;
    
    phi_km1km1 = -1/(3*hk^2)*((x_k - beta)^3 - (x_k - alpha)^3);
    phi_kkm1 = 1/hk^2*(-(beta^3 - alpha^3)/3 + ...
        (x_k + x_km1)*(beta^2 - alpha^2)/2 - x_km1*x_k*(beta - alpha));
    phi_kk = 1/(3*hk^2)*((beta - x_km1)^3 - (alpha - x_km1)^3);
    
    Ak = 0.5*[phi_km1km1, phi_kkm1; phi_kkm1, phi_kk];
      
    A_dt(kpos, kpos) = A_dt(kpos, kpos) - Ak;
    A_dt(kpos, M + kpos) = A_dt(kpos, M + kpos) + Ak;
    A_dt(M + kpos, kpos) = A_dt(M + kpos, kpos) - Ak;
    A_dt(M + kpos, M + kpos) = A_dt(M + kpos, M + kpos) + Ak;

end

