function A_git = GaInTedG1_AuxiliaryFunc(A_git, Kgs, x0, tg, nt, nx, ...
    GbL, GbR, Gb, hGb, gmin, gmax, w1, w2, pnltp)

global Ind kn tn tnm1

leKg = length(Kgs(1,:));
for k = 1:leKg                               
    
    x_km1 = Kgs(1, k);
    x_k = Kgs(2, k);
    hk = Kgs(3, k);
    
    x_kpos = find((x0 == x_k));
    kpos = (x_kpos - 1: x_kpos);
    
    alpha = max(gmin, x_km1);
    beta = min(gmax, x_k);
    alta = 0.5*(alpha + beta);  
    delx = beta - alpha;
    x_vec = [alpha alta beta];
        
    t_alpha = tg(k);
    t_beta = tg(k+1);
    t_alta = 0.5*(t_alpha + t_beta); 
    delt = abs(t_beta - t_alpha);
    t_vec = [t_alpha t_alta t_beta];
    
    delg = sqrt(delx^2 + delt^2); 
    
    w = zeros(1,3);    
    w(1) = delg/6;
    w(2) = 2*delg/3;
    w(3) = w(1);
        
    lamb1 = tn - t_vec;
    lamb2 = t_vec - tnm1;
       
    C11 = w.*lamb1.^2; 
    C12 = w.*lamb1.*lamb2;
    C22 = w.*lamb2.^2; 
    
    phi_km1 = x_k - x_vec;
    phi_k = x_vec - x_km1;
    phi_km1km1 = phi_km1.^2;
    phi_kkm1 = phi_k.*phi_km1;
    phi_kk = phi_k.^2; 
    
    Phi = 1/(kn^2)*[sum(C11) sum(C12); sum(C12) sum(C22)];
    Phi_km1 = 1/(hk*kn^2)*[C11*phi_km1' C12*phi_km1'; ...
                           C12*phi_km1' C22*phi_km1'];
    Phi_k = 1/(hk*kn^2)*[C11*phi_k' C12*phi_k'; C12*phi_k' C22*phi_k'];
    Phi_km1km1 = 1/(hk^2*kn^2)*[C11*phi_km1km1' C12*phi_km1km1'; ...
                                C12*phi_km1km1' C22*phi_km1km1'];
    Phi_kkm1 = 1/(hk^2*kn^2)*[C11*phi_kkm1' C12*phi_kkm1'; ...
                              C12*phi_kkm1' C22*phi_kkm1'];
    Phi_kk = 1/(hk^2*kn^2)*[C11*phi_kk' C12*phi_kk'; ...
                            C12*phi_kk' C22*phi_kk'];
  
    % The time term [u]v_d ------------------------------------------------
 
    % Recall that nt(in code) = -nt(in FE_form with n = n_1)
    if nt > 0       % Use v_1  
        
        A1k = nt*[Phi_km1km1 Phi_kkm1; Phi_kkm1 Phi_kk];       
        l = 1;
        for pos1 = kpos;
            m = 1;
            for pos2 = kpos;               
                A_git(Ind + pos1, Ind + pos2) = ...
                A_git(Ind + pos1, Ind + pos2) + A1k(l:l+1,m:m+1);
                m = m + 2;           
            end           
            l = l + 2;
        end
        
        A1Gbk = -nt*[Phi_km1; Phi_k];      
        l = 1;
        for pos1 = kpos;      
            A_git(Ind + pos1, Ind + Gb) = ...
            A_git(Ind + pos1, Ind + Gb) + A1Gbk(l:l+1,1:2);
            l = l + 2;
        end
    
    elseif nt < 0   % Use v_2
              
        A2k = nt*[Phi_km1 Phi_k];       
        m = 1;
        for pos2 = kpos;        
            A_git(Ind + Gb, Ind + pos2) = ...
            A_git(Ind + Gb, Ind + pos2) + A2k(1:2,m:m+1);
            m = m + 2;            
        end
        
        A2Gb = -nt*Phi;      
        A_git(Ind+Gb, Ind+Gb) = A_git(Ind+Gb, Ind+Gb) + A2Gb;
        
    end
    
    % The first space term <dn1u>[v] --------------------------------------   
    
    A1k = nx*w1/hk*[-Phi_km1, Phi_km1; -Phi_k, Phi_k];  
    l = 1;
    for pos1 = kpos;
        m = 1;
        for pos2 = kpos;      
            A_git(Ind + pos1, Ind + pos2) = ...
            A_git(Ind + pos1, Ind + pos2) + A1k(l:l+1,m:m+1);
            m = m + 2;
        end           
        l = l + 2;
    end
        
    A1GbL = -nx*w2/hGb*[Phi_km1; Phi_k];    
    l = 1;
    for pos1 = kpos;
        A_git(Ind + pos1, Ind + GbL) = ...
        A_git(Ind + pos1, Ind + GbL) + A1GbL(l:l+1,1:2);
        A_git(Ind + pos1, Ind + GbR) = ...
        A_git(Ind + pos1, Ind + GbR) - A1GbL(l:l+1,1:2);
        l = l + 2;
    end
    
    A2k = nx*w1/hk*[Phi, -Phi];   
    m = 1;
    for pos2 = kpos;
        A_git(Ind + Gb, Ind + pos2) = ...
        A_git(Ind + Gb, Ind + pos2) + A2k(1:2,m:m+1);
        m = m + 2;
    end
    
    A2GbL = nx*w2/hGb*Phi;   
    A_git(Ind + Gb, Ind + GbL) = A_git(Ind + Gb, Ind + GbL) + A2GbL;
    A_git(Ind + Gb, Ind + GbR) = A_git(Ind + Gb, Ind + GbR) - A2GbL;
    
    % The second space term, <dn1v>[u] ------------------------------------
    
    A1k = nx*w1/hk*[-Phi_km1, -Phi_k ; Phi_km1, Phi_k];   
    l = 1;
    for pos1 = kpos;
        m = 1;
        for pos2 = kpos;      
            A_git(Ind + pos1, Ind + pos2) = ...
            A_git(Ind + pos1, Ind + pos2) + A1k(l:l+1,m:m+1);
            m = m + 2;
        end           
        l = l + 2;
    end
    
    A1Gb = nx*w1/hk*[Phi; -Phi];
    l = 1;
    for pos1 = kpos;      
        A_git(Ind + pos1, Ind + Gb) = ...
        A_git(Ind + pos1, Ind + Gb) + A1Gb(l:l+1,1:2);
        l = l + 2;
    end
    
    A2k = nx*w2/hGb*[Phi_km1, Phi_k]; 
    m = 1;
    for pos2 = kpos;
        A_git(Ind + GbL, Ind + pos2) = ...
        A_git(Ind + GbL, Ind + pos2) - A2k(1:2,m:m+1);
        A_git(Ind + GbR, Ind + pos2) = ...
        A_git(Ind + GbR, Ind + pos2) + A2k(1:2,m:m+1);
        m = m + 2;
    end
    
    A2Gb = -nx*w2/hGb*Phi;
    A_git(Ind + GbL, Ind + Gb) = A_git(Ind + GbL, Ind + Gb) - A2Gb;
    A_git(Ind + GbR, Ind + Gb) = A_git(Ind + GbR, Ind + Gb) + A2Gb;

    % The stability term, [u][v] ------------------------------------------    
    
    hg = min(hk, hGb);
    %npen = abs(nx)*pnltp/hg;
    npen = pnltp/hg;
    
    A1k = npen*[Phi_km1km1, Phi_kkm1; Phi_kkm1, Phi_kk];
    l = 1;
    for pos1 = kpos;
        m = 1;
        for pos2 = kpos;      
            A_git(Ind + pos1, Ind + pos2) = ...
            A_git(Ind + pos1, Ind + pos2) + A1k(l:l+1,m:m+1);
            m = m + 2;
        end           
        l = l + 2;
    end
      
    A1Gb = -npen*[Phi_km1; Phi_k];  
    l = 1;
    for pos1 = kpos;      
        A_git(Ind + pos1, Ind + Gb) = ...
        A_git(Ind + pos1, Ind + Gb) + A1Gb(l:l+1,1:2);
        l = l + 2;
    end
    
    A2k = -npen*[Phi_km1, Phi_k];
    m = 1;
    for pos2 = kpos;
        A_git(Ind + Gb, Ind + pos2) = ...
        A_git(Ind + Gb, Ind + pos2) + A2k(1:2,m:m+1);
        m = m + 2;
    end
    
    A2Gb = npen*Phi;  
    A_git(Ind + Gb, Ind + Gb) = A_git(Ind + Gb, Ind + Gb) + A2Gb;
    
end