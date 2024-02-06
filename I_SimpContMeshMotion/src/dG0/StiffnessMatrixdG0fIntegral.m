function [A_sm, b_fint] = StiffnessMatrixdG0fIntegral(K0, KG, x0, G, ...
                                amin, amax, bmin, bmax)

global kn kn1ga kn1gb x0_init x0_fin I0 tnm1 tn h0 IG mutr

A_sm = zeros(I0 + IG);
Ak = zeros(2);

b_fint = zeros(I0 + IG, 1);
bk = zeros(2,1);

t_av = 0.5*(tn + tnm1);         % Time that is used in the quadrature

% Regular domains of U1 ---------------------------------------------------

% To the left of G
K1ra = K0(:, K0(2,:) < amin + h0);
leK1ra = length(K1ra(1,:));
for k = 1:leK1ra
    
    x_km1 = K1ra(1,k);
    x_k = K1ra(2,k);
    hk = K1ra(3, k);
    
    x_kpos = find((x0 == x_k));
    kpos = (x_kpos - 1: x_kpos);
    
    delx = min(amin, x_k) - x_km1;
    
    % A_sm
    Akel = kn/hk^2*delx;
    
    Ak = [Akel, -Akel ; -Akel, Akel];
    
    A_sm(kpos, kpos) = A_sm(kpos, kpos) + Ak;
    
    % b_fint
    alpha_kn = 0.5*(f_func(x_km1,t_av) + f_func(x_k,t_av));
    beta_kn = 0.5*(f_func(x_km1,t_av)*x_km1 + f_func(x_k,t_av)*x_k);
    
    bk = 1/hk*delx*kn*[x_k*alpha_kn - beta_kn; beta_kn - x_km1*alpha_kn];
    
    b_fint(kpos) = b_fint(kpos) + bk;
    
end

% To the right of G
K1rb = K0(:, bmax - h0 < K0(1,:));
leK1rb = length(K1rb(1,:));
for k = 1:leK1rb
    
    x_km1 = K1rb(1,k);
    x_k = K1rb(2,k);
    hk = K1rb(3, k);
    
    x_kpos = find((x0 == x_k));
    kpos = (x_kpos - 1: x_kpos);
    
    delx = x_k - max(bmax, x_km1);
    
    % A_sm
    Akel = kn/hk^2*delx;
    
    Ak = [Akel, -Akel ; -Akel, Akel];
    
    A_sm(kpos, kpos) = A_sm(kpos, kpos) + Ak;
    
    % b_fint
    alpha_kn = 0.5*(f_func(x_km1,t_av) + f_func(x_k,t_av));
    beta_kn = 0.5*(f_func(x_km1,t_av)*x_km1 + f_func(x_k,t_av)*x_k);
    
    bk = 1/hk*delx*kn*[x_k*alpha_kn - beta_kn; beta_kn - x_km1*alpha_kn];
    
    b_fint(kpos) = b_fint(kpos) + bk;
    
end

% Gamma domains of U1 -----------------------------------------------------

% To the left of G
K1ga = K0(:, amin - h0 < K0(1,:) & K0(2,:) < amax + h0);
leK1ga = length(K1ga(1,:));
for k = 1:leK1ga
    
    x_km1 = K1ga(1,k);
    x_k = K1ga(2,k);
    hk = K1ga(3, k);
    
    x_kpos = find((x0 == x_k));
    kpos = (x_kpos - 1: x_kpos);
    
    delx = min(amax, x_k) - max(amin, x_km1);
    intdom = 0.5*delx*(kn1ga(1,k) + kn1ga(1,k+1));
    
    % A_sm
    Akel = intdom/hk^2; 
    
    Ak = [Akel, -Akel ; -Akel, Akel];
    
    A_sm(kpos, kpos) = A_sm(kpos, kpos) + Ak;
    
    % b_fint
    t_km1q = kn1ga(2,k);
    t_kq = kn1ga(2,k+1);
    
    alpha_kn = 0.5*(f_func(x_km1,t_km1q) + f_func(x_k,t_kq));
    beta_kn = 0.5*(f_func(x_km1,t_km1q)*x_km1 + f_func(x_k,t_kq)*x_k);
    
    bk = intdom/hk*[x_k*alpha_kn - beta_kn; beta_kn - x_km1*alpha_kn];
    
    b_fint(kpos) = b_fint(kpos) + bk;
    
end

% To the right of G
K1gb = K0(:, bmin - h0 < K0(1,:) & K0(2,:) < bmax + h0);
leK1gb = length(K1gb(1,:));
for k = 1:leK1gb
    
    x_km1 = K1gb(1,k);
    x_k = K1gb(2,k);
    hk = K1gb(3, k);
    
    x_kpos = find((x0 == x_k));
    kpos = (x_kpos - 1: x_kpos);
    
    delx = min(bmax, x_k) - max(bmin, x_km1);
    intdom = 0.5*delx*(kn1gb(1,k) + kn1gb(1,k+1));
    
    % A_sm
    Akel = intdom/hk^2;  
    
    Ak = [Akel, -Akel ; -Akel, Akel];
    
    A_sm(kpos, kpos) = A_sm(kpos, kpos) + Ak;    
    
    % b_fint
    t_km1q = kn1gb(2,k);
    t_kq = kn1gb(2,k+1);
    
    alpha_kn = 0.5*(f_func(x_km1,t_km1q) + f_func(x_k,t_kq));
    beta_kn = 0.5*(f_func(x_km1,t_km1q)*x_km1 + f_func(x_k,t_kq)*x_k);
    
    bk = intdom/hk*[x_k*alpha_kn - beta_kn; beta_kn - x_km1*alpha_kn];
    
    b_fint(kpos) = b_fint(kpos) + bk;
end

% Domains of U2 -----------------------------------------------------------

K2r = KG(:, x0_init <= KG(1,:) & KG(2,:) <= x0_fin);
leK2r = length(K2r(1,:));
for k = 1:leK2r
    
    x_km1 = K2r(1,k);
    x_k = K2r(2,k);
    hk = K2r(3, k);
    
    x_kpos = find((G == x_k));
    kpos = I0 + (x_kpos - 1: x_kpos);
    
    % A_sm
    Akel = kn/hk;
    
    Ak = [Akel, -Akel ; -Akel, Akel];
    
    A_sm(kpos, kpos) = A_sm(kpos, kpos) + Ak;
    
    % b_fint
    x_km1q = x_km1 - mutr*(tn - t_av);
    x_kq = x_k - mutr*(tn - t_av);
    
    alpha_kn = 0.5*(f_func(x_km1q,t_av) + f_func(x_kq,t_av));
    beta_kn = 0.5*(f_func(x_km1q,t_av)*x_km1q + f_func(x_kq,t_av)*x_kq);
    
    bk = kn*[x_kq*alpha_kn - beta_kn; beta_kn - x_km1q*alpha_kn];
    
    b_fint(kpos) = b_fint(kpos) + bk;
    
end

% figure(1)
% test = K1ra(1:2,:);
% plot(test,tn*ones(2, length(test(1,:))),'go--');
% hold on
% test = K1ga(1:2,:);
% plot(test,tn*ones(2, length(test(1,:))),'go--');

    
    
    
    
    
    
    
    
    
    