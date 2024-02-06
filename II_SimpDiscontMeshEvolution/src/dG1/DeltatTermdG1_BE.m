function A_dt = DeltatTermdG1_BE(K0, KG, x0, G)

global x0_init x0_fin I0 leA M 

A_dt = zeros(leA);

xGl = G(1);
xGr = G(end);

% Domains of U1 ---------------------------------------------------

% To the left of G
K1ra = K0(:, K0(1,:) < xGl);
leK1ra = length(K1ra(1,:));
for k = 1:leK1ra
    
    x_km1 = K1ra(1,k);
    x_k = K1ra(2,k);
    hk = K1ra(3, k);
    
    x_kpos = find((x0 == x_k));
    kpos = (x_kpos - 1: x_kpos);
    
    alpha = x_km1;
    beta = min(x_k, xGl);
    
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
K1rb = K0(:, xGr < K0(2,:)); 
leK1rb = length(K1rb(1,:));
for k = 1:leK1rb
    
    x_km1 = K1rb(1,k);
    x_k = K1rb(2,k);
    hk = K1rb(3, k);
    
    x_kpos = find((x0 == x_k));
    kpos = (x_kpos - 1: x_kpos);
    
    alpha = max(x_km1, xGr);
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

