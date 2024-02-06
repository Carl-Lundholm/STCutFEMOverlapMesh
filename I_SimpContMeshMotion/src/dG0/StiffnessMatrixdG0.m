function A_sm = StiffnessMatrixdG0(K0, KG, x0, G, ...
                                amin, amax, bmin, bmax)

global kn kn1ga kn1gb x0_init x0_fin I0 IG 

A_sm = zeros(I0 + IG);

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
    
    delx = min(amin, x_k) - x_km1;

    Akel = kn/hk^2*delx;
    
    Ak = [Akel, -Akel ; -Akel, Akel];
    
    A_sm(kpos, kpos) = A_sm(kpos, kpos) + Ak;
   
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
    
    delx = x_k - max(bmax, x_km1);

    Akel = kn/hk^2*delx;
    
    Ak = [Akel, -Akel ; -Akel, Akel];
    
    A_sm(kpos, kpos) = A_sm(kpos, kpos) + Ak;
    
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
    
    delx = min(amax, x_k) - max(amin, x_km1);
    intdom = 0.5*delx*(kn1ga(1,k) + kn1ga(1,k+1));

    Akel = intdom/hk^2; 
    
    Ak = [Akel, -Akel ; -Akel, Akel];
    
    A_sm(kpos, kpos) = A_sm(kpos, kpos) + Ak;
    
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
    
    delx = min(bmax, x_k) - max(bmin, x_km1);
    intdom = 0.5*delx*(kn1gb(1,k) + kn1gb(1,k+1));

    Akel = intdom/hk^2;  
    
    Ak = [Akel, -Akel ; -Akel, Akel];
    
    A_sm(kpos, kpos) = A_sm(kpos, kpos) + Ak;    

end

% Domains of U2 -----------------------------------------------------------

K2r = KG(:, x0_init <= KG(1,:) & KG(2,:) <= x0_fin);
leK2r = length(K2r(1,:));
for k = 1:leK2r
    
    x_k = K2r(2,k);
    hk = K2r(3, k);
    
    x_kpos = find((G == x_k));
    kpos = I0 + (x_kpos - 1: x_kpos);
    
    Akel = kn/hk;
    
    Ak = [Akel, -Akel ; -Akel, Akel];
    
    A_sm(kpos, kpos) = A_sm(kpos, kpos) + Ak;

end

% figure(1)
% test = K1ra(1:2,:);
% plot(test,tn*ones(2, length(test(1,:))),'go--');
% hold on
% test = K1ga(1:2,:);
% plot(test,tn*ones(2, length(test(1,:))),'go--');

    
    
    
    
    
    
    
    
    
    