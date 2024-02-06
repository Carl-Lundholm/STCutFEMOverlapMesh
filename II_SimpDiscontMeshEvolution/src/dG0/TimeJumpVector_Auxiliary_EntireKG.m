function b_k = TimeJumpVector_Auxiliary_EntireKG(Kp, Unm1, K0s, KGs, x0, G)

global I0

% Checks which simplices in mesh_composition at tnm1_minus that entire 
% (uncut) KG-simplex Kp intersects and computes the local contribution to 
% rhs-vector.

% Simplex Kp belongs to mesh composition at time tnm1_plus
x_km1 = Kp(1);
x_k = Kp(2); 

% The left and right endpoints xGl and xGr belong to G(tnm1_minus)
xGl = G(1);
xGr = G(end);

% Check intersection and compute the contribution
if x_k <= xGl || x_km1 >= xGr % Kp only intersects K0-simplices
    
    b_k = [0; 0];
   
    % K0sCC, CC = Covered and Cut by Kp
    K0sCC = K0s(:, x_km1 < K0s(2,:) & K0s(1,:) < x_k);   
    for Km = K0sCC
        
        x_lm1 = Km(1);
        x_l = Km(2);
        x_lpos = find((x0 == x_l));
        lpos = (x_lpos - 1: x_lpos);
        
        x_le = max(x_lm1, x_km1);   % Modify left endpoint if needed
        x_re = min(x_l, x_k);       % Modify right endpoint if needed
        
        b_kloc = TimeJumpVector_Auxiliary_LocalofLocal_DiffMesh(x_le, ...
                    x_re, Unm1(lpos(1)), Unm1(lpos(2)), Kp, Km);
        
        b_k = b_k + b_kloc;
        
    end
    
elseif x_km1 < xGl && xGl < x_k % Cut: Kp contains xGl(tnm1_minus)
    
    b_k = [0; 0];
    
    % Omega_1-part
    % K0sCC, CC = Covered and Cut by Kp
    K0sCC = K0s(:, x_km1 < K0s(2,:) & K0s(1,:) < xGl);   
    for Km = K0sCC
        
        x_lm1 = Km(1);
        x_l = Km(2);
        x_lpos = find((x0 == x_l));
        lpos = (x_lpos - 1: x_lpos);
        
        x_le = max(x_lm1, x_km1);   % Modify left endpoint if needed
        x_re = min(x_l, xGl);       % Modify right endpoint if needed
        
        b_kloc = TimeJumpVector_Auxiliary_LocalofLocal_DiffMesh(x_le, ...
                    x_re, Unm1(lpos(1)), Unm1(lpos(2)), Kp, Km);
        
        b_k = b_k + b_kloc;
        
    end
    
    % Omega_2/G-part
    % KGsCC, CC = Covered and Cut by Kp
    KGsCC = KGs(:, KGs(1,:) < x_k);
    for Km = KGsCC
        
        x_lm1 = Km(1);
        x_l = Km(2);
        x_lpos = find((G == x_l));
        lpos = I0 + (x_lpos - 1: x_lpos);
        
        x_le = x_lm1;
        x_re = min(x_l, x_k);       % Modify right endpoint if needed
        
        b_kloc = TimeJumpVector_Auxiliary_LocalofLocal_DiffMesh(x_le, ...
                    x_re, Unm1(lpos(1)), Unm1(lpos(2)), Kp, Km);
        
        b_k = b_k + b_kloc;
        
    end
        
elseif x_km1 < xGr && xGr < x_k % Cut: Kp contains xGr
    
    b_k = [0; 0];
    
    % Omega_2/G-part
    % KGsCC, CC = Covered and Cut by Kp
    KGsCC = KGs(:, KGs(2,:) > x_km1);
    for Km = KGsCC
        
        x_lm1 = Km(1);
        x_l = Km(2);
        x_lpos = find((G == x_l));
        lpos = I0 + (x_lpos - 1: x_lpos);
        
        x_le = max(x_lm1, x_km1); % Modify left endpoint if needed
        x_re = x_l;
       
        b_kloc = TimeJumpVector_Auxiliary_LocalofLocal_DiffMesh(x_le, ...
                    x_re, Unm1(lpos(1)), Unm1(lpos(2)), Kp, Km);
        
        b_k = b_k + b_kloc;
        
    end
    
    % Omega_1-part
    % K0sCC, CC = Covered and Cut by Kp
    K0sCC = K0s(:, xGr < K0s(2,:) & K0s(1,:) < x_k);
    for Km = K0sCC
        
        x_lm1 = Km(1);
        x_l = Km(2);
        x_lpos = find((x0 == x_l));
        lpos = (x_lpos - 1: x_lpos);
        
        x_le = max(x_lm1, xGr); % Modify left endpoint if needed
        x_re = min(x_l, x_k);
       
        b_kloc = TimeJumpVector_Auxiliary_LocalofLocal_DiffMesh(x_le, ...
                    x_re, Unm1(lpos(1)), Unm1(lpos(2)), Kp, Km);
        
        b_k = b_k + b_kloc;
        
    end
      

else % Covered: Kp lies completely in Gnm1_minus (xGl<=x_km1 & x_k<=xGr)
   
    b_k = [0; 0];
   
    % KGsCC, CC = Covered and Cut by Kp
    KGsCC = KGs(:, x_km1 < KGs(2,:) & KGs(1,:) < x_k);   
    for Km = KGsCC
        
        x_lm1 = Km(1);
        x_l = Km(2);
        x_lpos = find((G == x_l));
        lpos = I0 + (x_lpos - 1: x_lpos);
        
        x_le = max(x_lm1, x_km1);   % Modify left endpoint if needed
        x_re = min(x_l, x_k);       % Modify right endpoint if needed
        
        b_kloc = TimeJumpVector_Auxiliary_LocalofLocal_DiffMesh(x_le, ...
                    x_re, Unm1(lpos(1)), Unm1(lpos(2)), Kp, Km);
        
        b_k = b_k + b_kloc;
        
    end
    
end



