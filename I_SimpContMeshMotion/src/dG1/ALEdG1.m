function A_ale = ALEdG1(KG, G, mutr, kn)

global I0 leA M x0_init x0_fin 

A_ale = zeros(leA);

K2r = KG(:, x0_init <= KG(1,:) & KG(2,:) <= x0_fin);
leK2r = length(K2r(1,:));

for k = 1:leK2r
    
    x_k = K2r(2,k);
    
    x_kpos = find((G == x_k));
    kpos = I0 + (x_kpos - 1: x_kpos);
    
    Ak = mutr*kn/(12)*[1 -1; 1 -1];
      
    A_ale(kpos, kpos) = A_ale(kpos, kpos) + 2*Ak;
    A_ale(kpos, M + kpos) = A_ale(kpos, M + kpos) + Ak;
    A_ale(M + kpos, kpos) = A_ale(M + kpos, kpos) + Ak;
    A_ale(M + kpos, M + kpos) = A_ale(M + kpos, M + kpos) + 2*Ak;

end
