function A_sm = StiffnessMatrixdG1_BE(K0, KG, x0, G)

global kn x0_init x0_fin I0 leA M

A_sm = zeros(leA);

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
    
    Ak = (beta - alpha)*kn/(6*hk^2)*[1 -1; -1 1];
      
    A_sm(kpos, kpos) = A_sm(kpos, kpos) + 2*Ak;
    A_sm(kpos, M + kpos) = A_sm(kpos, M + kpos) + Ak;
    A_sm(M + kpos, kpos) = A_sm(M + kpos, kpos) + Ak;
    A_sm(M + kpos, M + kpos) = A_sm(M + kpos, M + kpos) + 2*Ak;
   
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
    
    Ak = (beta - alpha)*kn/(6*hk^2)*[1 -1; -1 1];
      
    A_sm(kpos, kpos) = A_sm(kpos, kpos) + 2*Ak;
    A_sm(kpos, M + kpos) = A_sm(kpos, M + kpos) + Ak;
    A_sm(M + kpos, kpos) = A_sm(M + kpos, kpos) + Ak;
    A_sm(M + kpos, M + kpos) = A_sm(M + kpos, M + kpos) + 2*Ak;
    
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

    
    
    
    
    
    
    
    
    
    