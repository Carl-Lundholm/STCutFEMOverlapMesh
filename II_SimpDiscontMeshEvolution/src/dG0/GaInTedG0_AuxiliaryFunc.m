function A_git = GaInTedG0_AuxiliaryFunc(A_git, Kgs, x0, C, tg, nt, nx, ...
    GbL, GbR, Gb, hGb, gmin, gmax, w1, w2, pnltp)

leKgs = length(Kgs(1,:));

for k = 1:leKgs                              
                                    
    x_km1 = Kgs(1, k);
    x_k = Kgs(2, k);
    hk = Kgs(3, k);
    
    t_k = tg(k + 1);
    t_km1 = tg(k);
    
    x_kpos = find((x0 == x_k));
    kpos = (x_kpos - 1: x_kpos);
    
    alpha0 = max(gmin, x_km1);
    beta0 = min(gmax, x_k);
    
    alpha = alpha0 - x_km1;
    beta = beta0 - x_km1;
    
    delx = beta0 - alpha0;
    delt = t_k - t_km1;
    
    gamma_k = abs(delt*C);
 
    I_km1 = gamma_k*(1 - (delx/2 + alpha)/hk);
    I_k = gamma_k/hk*(delx/2 + alpha);
    I_km1km1 = gamma_k*(1 - (alpha + beta)/hk + ...
                           (delx^2/3 + alpha*beta)/hk^2);
    I_kk = gamma_k/(3*hk^2)*(alpha^2 + alpha*beta + beta^2);
    I_kkm1 = I_k - I_kk;
    
    % The time term [u]v_d
 
    if nt > 0       % Use v_1
    
        A1k = nt*[I_km1km1, I_kkm1 ; I_kkm1, I_kk];
        A1Gb = -nt*[I_km1; I_k];
        
        A_git(kpos, kpos) = A_git(kpos, kpos) + A1k;
        A_git(kpos, Gb) = A_git(kpos, Gb) + A1Gb;
    
    elseif nt < 0   % Use v_2
        
        A2k = nt*[I_km1, I_k];
        A2Gb = -nt*gamma_k;
        
        A_git(Gb, kpos) = A_git(Gb, kpos) + A2k;
        A_git(Gb, Gb) = A_git(Gb, Gb) + A2Gb;
        
    end
    
    % The first space term <dn1u>[v]   
    
    A1k = nx*w1/hk*[-I_km1, I_km1 ; -I_k, I_k];
    A1GbL = -nx*w2/hGb*[I_km1; I_k];
    A1GbR = -A1GbL;
    
    A_git(kpos, kpos) = A_git(kpos, kpos) + A1k;
    A_git(kpos, GbL) = A_git(kpos, GbL) + A1GbL;
    A_git(kpos, GbR) = A_git(kpos, GbR) + A1GbR;
    
    A2k = nx*w1*gamma_k/hk*[1, -1]; 
    A2GbL = nx*w2/hGb*gamma_k;
    A2GbR = -A2GbL;
    
    A_git(Gb, kpos) =  A_git(Gb, kpos) + A2k;
    A_git(Gb, GbL) = A_git(Gb, GbL) + A2GbL;
    A_git(Gb, GbR) = A_git(Gb, GbR) + A2GbR;
    
    % The second space term, <dn1v>[u]
    
    A1k = nx*w1/hk*[-I_km1, -I_k ; I_km1, I_k];
    A1Gb = nx*w1*gamma_k/hk*[1; -1];
    
    A_git(kpos, kpos) = A_git(kpos, kpos) + A1k;
    A_git(kpos, Gb) = A_git(kpos, Gb) + A1Gb;
    
    A2k = nx*w2/hGb*[I_km1, I_k];  
    A2Gb = -nx*w2/hGb*gamma_k;
    
    A_git(GbL, kpos) = A_git(GbL, kpos) - A2k;
    A_git(GbR, kpos) = A_git(GbR, kpos) + A2k;
    A_git(GbL, Gb) = A_git(GbL, Gb) - A2Gb;
    A_git(GbR, Gb) = A_git(GbR, Gb) + A2Gb;

    % The stability term, [u][v]    
    
    hg = min(hk, hGb);
    npen = pnltp/hg;
    
    A1k = npen*[I_km1km1, I_kkm1; I_kkm1, I_kk];
    A1Gb = -npen*[I_km1; I_k];  
    
    A_git(kpos, kpos) = A_git(kpos, kpos) + A1k;
    A_git(kpos, Gb) = A_git(kpos, Gb) + A1Gb;
    
    A2k = -npen*[I_km1, I_k];
    A2Gb = npen*gamma_k;  
    
    A_git(Gb, kpos) = A_git(Gb, kpos) + A2k;
    A_git(Gb, Gb) = A_git(Gb, Gb) + A2Gb;
    
end