function A_pnt = PeNaTedG1_AuxiliaryFunc(A_pnt, K0, x0, n1, nx, tg, ...
    GbL, GbR, hGb, gmin, gmax, pnltp_nab)

global tnm1 tn Ind kn

Kg = K0(:, K0(1,:) < gmax & gmin < K0(2,:));
leKg = length(Kg(1,:));

for k = 1:leKg                              
                                   
    x_km1 = Kg(1, k);
    x_k = Kg(2, k);
    hk = Kg(3, k);
    
    x_kpos = find((x0 == x_k));
    kpos = (x_kpos - 1: x_kpos);
    
    alpha = max(gmin, x_km1);
    beta = min(gmax, x_k);
    alta = 0.5*(alpha + beta);
    x_vec = [alpha alta beta];
    
    t_alpha = tg(k);
    t_beta = tg(k + 1);
    t_alta = 0.5*(t_alpha + t_beta);
    t_vec = [t_alpha t_alta t_beta];
    
    deltg = abs(t_beta - t_alpha); 
    
    w = zeros(1,3);    
    w(1) = deltg/6;
    w(2) = 2*deltg/3;
    w(3) = w(1);
        
    lamb1 = tn - t_vec;
    lamb2 = t_vec - tnm1;
       
    C11 = w.*lamb1.^2; 
    C12 = w.*lamb1.*lamb2;
    C22 = w.*lamb2.^2; 
    
    if n1 > 0       % Overlap regions to the left of G 
        delOmegaO = x_k - x_vec;        
    
    else            % Overlap regions to the right of G     
        delOmegaO = x_vec - x_km1;
        
    end
    
    Phi = [C11*delOmegaO' C12*delOmegaO'; C12*delOmegaO' C22*delOmegaO'];
       
    % The derivative stability term, [du][dv]
    
    fac = pnltp_nab/kn^2; 
    
    A1k = fac/hk^2*[Phi, -Phi ; -Phi, Phi];
    A1GbL = fac/(hk*hGb)*[-Phi; Phi];    
    l = 1;
    for pos1 = kpos;       
        m = 1;
        for pos2 = kpos;
            A_pnt(Ind + pos1, Ind + pos2) = ...
            A_pnt(Ind + pos1, Ind + pos2) + A1k(l:l+1, m:m+1);
            m = m + 2;
        end
        A_pnt(Ind + pos1, Ind + GbL) = ...
        A_pnt(Ind + pos1, Ind + GbL) + A1GbL(l:l+1, 1:2);
        A_pnt(Ind + pos1, Ind + GbR) = ...
        A_pnt(Ind + pos1, Ind + GbR) - A1GbL(l:l+1, 1:2);
        l = l + 2;
    end 
        
    A2k = fac/(hGb*hk)*[Phi, -Phi];
    m = 1;
    for pos2 = kpos;      
        A_pnt(Ind + GbL, Ind + pos2) = ...
        A_pnt(Ind + GbL, Ind + pos2) - A2k(1:2, m:m+1);
    
        A_pnt(Ind + GbR, Ind + pos2) = ...
        A_pnt(Ind + GbR, Ind + pos2) + A2k(1:2, m:m+1);
    
        m = m + 2;      
    end
       
%     A2Gb = fac/(hGb^2)*Phi;    
%     A_pnt(GbL, GbL) = A_pnt(GbL, GbL) + A2Gb;
%     A_pnt(GbR, GbL) = A_pnt(GbR, GbL) - A2Gb;
%     A_pnt(GbL, GbR) = A_pnt(GbL, GbR) - A2Gb;
%     A_pnt(GbR, GbR) = A_pnt(GbR, GbR) + A2Gb;    
    A2Gb_mat = fac/(hGb^2)*[Phi, -Phi; -Phi, Phi];
    l = 1;
    for pos1 = [GbL,GbR]
        m = 1;
        for pos2 = [GbL,GbR]
            A_pnt(Ind + pos1, Ind + pos2) = ...
            A_pnt(Ind + pos1, Ind + pos2) + A2Gb_mat(l:l+1, m:m+1);
            m = m + 2;
        end
        l = l + 2;
    end
    
end