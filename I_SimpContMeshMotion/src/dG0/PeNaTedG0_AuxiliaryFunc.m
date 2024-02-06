function A_pnt = PeNaTedG0_AuxiliaryFunc(A_pnt, K0, x0, n1, nx, tg, ...
    GbL, GbR, hGb, gmin, gmax, pnltp_nab)

Kg = K0(:, K0(1,:) < gmax & gmin < K0(2,:));
leKg = length(Kg(1,:));

for k = 1:leKg                              
                                   
    x_km1 = Kg(1, k);
    x_k = Kg(2, k);
    hk = Kg(3, k);
    
    t_k = tg(k + 1);
    t_km1 = tg(k);
    
    x_kpos = find((x0 == x_k));
    kpos = (x_kpos - 1: x_kpos);
    
    alpha0 = max(gmin, x_km1);
    beta0 = min(gmax, x_k);
    
    if n1 > 0       % Overlap region to the left of G
        delOmegaO = 2*x_k - alpha0 - beta0;        
    
    else            % Overlap region to the right of G     
        delOmegaO = alpha0 + beta0 - 2*x_km1;
        
    end
    
%     delOmegaO = hk*2;
    
    delt = abs(t_k - t_km1);
%     delt = 1;
    
    intdom = 0.5*delOmegaO*delt;
       
    % The derivative stability term, [du][dv]
    
    %fac = abs(nx)*pnltp_nab*intdom;
    fac = pnltp_nab*intdom; 
    
    A1k = fac/hk^2*[1, -1 ; -1, 1];
    A1GbL = fac/(hk*hGb)*[-1; 1];
    A1GbR = -A1GbL;
    
    A_pnt(kpos, kpos) = A_pnt(kpos, kpos) + A1k;
    A_pnt(kpos, GbL) = A_pnt(kpos, GbL) + A1GbL;
    A_pnt(kpos, GbR) = A_pnt(kpos, GbR) + A1GbR;
    
    A2k = fac/(hGb*hk)*[1, -1];
    A2Gb = fac/hGb^2;
    A2Gb_mat = A2Gb*[1, -1; -1, 1];
    
    A_pnt(GbL, kpos) = A_pnt(GbL, kpos) - A2k;
    A_pnt(GbR, kpos) = A_pnt(GbR, kpos) + A2k;
%     A_pnt(GbL, GbL) = A_pnt(GbL, GbL) + A2Gb;
%     A_pnt(GbR, GbL) = A_pnt(GbR, GbL) - A2Gb;
%     A_pnt(GbL, GbR) = A_pnt(GbL, GbR) - A2Gb;
%     A_pnt(GbR, GbR) = A_pnt(GbR, GbR) + A2Gb;
    
    A_pnt(GbL:GbR, GbL:GbR) = A_pnt(GbL:GbR, GbL:GbR) + A2Gb_mat;
    
    
end