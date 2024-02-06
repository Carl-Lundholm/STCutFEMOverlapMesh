function A_sm = StiffnessMatrixdG1(K0, KG, x0, G, ...
                                amin, amax, bmin, bmax, t1ga, t1gb)

global kn x0_init x0_fin I0 leA M tnm1 tn

A_sm = zeros(leA);

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
    
    Ak = (beta - alpha)*kn/(6*hk^2)*[1 -1; -1 1];
      
    A_sm(kpos, kpos) = A_sm(kpos, kpos) + 2*Ak;
    A_sm(kpos, M + kpos) = A_sm(kpos, M + kpos) + Ak;
    A_sm(M + kpos, kpos) = A_sm(M + kpos, kpos) + Ak;
    A_sm(M + kpos, M + kpos) = A_sm(M + kpos, M + kpos) + 2*Ak;
   
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
    
    Ak = (beta - alpha)*kn/(6*hk^2)*[1 -1; -1 1];
      
    A_sm(kpos, kpos) = A_sm(kpos, kpos) + 2*Ak;
    A_sm(kpos, M + kpos) = A_sm(kpos, M + kpos) + Ak;
    A_sm(M + kpos, kpos) = A_sm(M + kpos, kpos) + Ak;
    A_sm(M + kpos, M + kpos) = A_sm(M + kpos, M + kpos) + 2*Ak;
    
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
    
    phi_mat = (beta - alpha)/(hk^2)*[1 -1; -1 1];
    
    beta = 0.5*(alpha + beta);
    
    phi_matgm = (beta - alpha)/(hk^2)*[1 -1; -1 1];
    
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
    
    C11gm = wgm*lamb1gm*lamb1gm;
    C11 = (wgf + wf)*lamb1f*lamb1f + wm*lamb1m*lamb1m + we*lamb1e*lamb1e; 
    
    C12gm = wgm*lamb1gm*lamb2gm;
    C12 = (wgf + wf)*lamb1f*lamb2f + wm*lamb1m*lamb2m + we*lamb1e*lamb2e;
    
    C22gm = wgm*lamb2gm*lamb2gm;
    C22 = (wgf + wf)*lamb2f*lamb2f + wm*lamb2m*lamb2m + we*lamb2e*lamb2e;
    
    Ak11 = (kn^-2)*(C11gm*phi_matgm + C11*phi_mat);
    Ak12 = (kn^-2)*(C12gm*phi_matgm + C12*phi_mat);
    Ak22 = (kn^-2)*(C22gm*phi_matgm + C22*phi_mat);
      
    A_sm(kpos, kpos) = A_sm(kpos, kpos) + Ak11;
    A_sm(kpos, M + kpos) = A_sm(kpos, M + kpos) + Ak12;
    A_sm(M + kpos, kpos) = A_sm(M + kpos, kpos) + Ak12;
    A_sm(M + kpos, M + kpos) = A_sm(M + kpos, M + kpos) + Ak22;
    
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
    
    phi_mat = (beta - alpha)/(hk^2)*[1 -1; -1 1];
    
    alpha = 0.5*(alpha + beta);
    
    phi_matgm = (beta - alpha)/(hk^2)*[1 -1; -1 1];
    
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
    
    C11gm = wgm*lamb1gm*lamb1gm;
    C11 = (wgf + wf)*lamb1f*lamb1f + wm*lamb1m*lamb1m + we*lamb1e*lamb1e; 
    
    C12gm = wgm*lamb1gm*lamb2gm;
    C12 = (wgf + wf)*lamb1f*lamb2f + wm*lamb1m*lamb2m + we*lamb1e*lamb2e;
    
    C22gm = wgm*lamb2gm*lamb2gm;
    C22 = (wgf + wf)*lamb2f*lamb2f + wm*lamb2m*lamb2m + we*lamb2e*lamb2e;
    
    Ak11 = (kn^-2)*(C11gm*phi_matgm + C11*phi_mat);
    Ak12 = (kn^-2)*(C12gm*phi_matgm + C12*phi_mat);
    Ak22 = (kn^-2)*(C22gm*phi_matgm + C22*phi_mat);
      
    A_sm(kpos, kpos) = A_sm(kpos, kpos) + Ak11;
    A_sm(kpos, M + kpos) = A_sm(kpos, M + kpos) + Ak12;
    A_sm(M + kpos, kpos) = A_sm(M + kpos, kpos) + Ak12;
    A_sm(M + kpos, M + kpos) = A_sm(M + kpos, M + kpos) + Ak22;  

end

% Domains of U2 -----------------------------------------------------------

K2r = KG(:, x0_init <= KG(1,:) & KG(2,:) <= x0_fin);
leK2r = length(K2r(1,:));
for k = 1:leK2r
    
    x_km1 = K2r(1,k);
    x_k = K2r(2,k);
    hk = K2r(3,k);
    
    x_kpos = find((G == x_k));
    kpos = I0 + (x_kpos - 1: x_kpos);
    
    alpha = x_km1;
    beta = x_k;
    
    Ak = (beta - alpha)*kn/(6*hk^2)*[1 -1; -1 1];
      
    A_sm(kpos, kpos) = A_sm(kpos, kpos) + 2*Ak;
    A_sm(kpos, M + kpos) = A_sm(kpos, M + kpos) + Ak;
    A_sm(M + kpos, kpos) = A_sm(M + kpos, kpos) + Ak;
    A_sm(M + kpos, M + kpos) = A_sm(M + kpos, M + kpos) + 2*Ak;

end

% figure(1)
% test = K1ra(1:2,:);
% plot(test,tn*ones(2, length(test(1,:))),'go--');
% hold on
% test = K1ga(1:2,:);
% plot(test,tn*ones(2, length(test(1,:))),'go--');

    
    
    
    
    
    
    
    
    
    