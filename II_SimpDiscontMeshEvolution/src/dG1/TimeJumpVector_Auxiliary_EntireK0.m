function b_k = TimeJumpVector_Auxiliary_EntireK0(Kp, kpos, Unm1, KGs, G)

global I0

% Checks which simplices in mesh_composition at tnm1_minus that entire 
% (uncut) K0-simplex Kp intersects and computes the local contribution to 
% rhs-vector.

% Note that this is m-file is a special case and somewhat cheaper 
% version of the m-file "TimeJumpVector_Auxiliary_CutK0".

% Simplex Kp belongs to mesh composition at time tnm1_plus
x_km1 = Kp(1);
x_k = Kp(2);
h_k = Kp(3); 

% The left and right endpoints xGl and xGr belong to G(tnm1_minus)
xGl = G(1);
xGr = G(end);

% Check intersection and compute the contribution
if x_k <= xGl || x_km1 >= xGr % Uncovered: Kp is unchanged (Standard case?)
    
    U_km1 = Unm1(kpos(1));
    U_k = Unm1(kpos(2));
    
    b_k = h_k/6*[2*U_km1 + U_k; U_km1 + 2*U_k];
    
elseif x_km1 < xGl && xGl < x_k % Cut: Kp contains xGl(tnm1_minus)
    
    b_k = [0; 0];
    
    % Omega_1-part
    b_kloc = TimeJumpVector_Auxiliary_LocalofLocal_SameMesh(x_km1, xGl, ...
                    Unm1(kpos(1)), Unm1(kpos(2)), Kp);
                
    b_k = b_k + b_kloc;    
    
    % Omega_2/G-part
    % KGsCC, CC = Covered and Cut by Kp
    KGsCC = KGs(:, KGs(1,:) < x_k);
    for Km = KGsCC
        
        x_lm1 = Km(1);
        x_l = Km(2);
        x_lpos = find((G == x_l));
        lpos = I0 + (x_lpos - 1: x_lpos);
        
        x_le = x_lm1;
        x_re = min(x_l, x_k); % Modify right endpoint if needed
        
        b_kloc = TimeJumpVector_Auxiliary_LocalofLocal_DiffMesh(x_le, ...
                    x_re, Unm1(lpos(1)), Unm1(lpos(2)), Kp, Km);
        
        b_k = b_k + b_kloc;
    end
        
elseif x_km1 < xGr && xGr < x_k % Cut: Kp contains xG_right
    
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
    b_kloc = TimeJumpVector_Auxiliary_LocalofLocal_SameMesh(xGr, x_k, ...
                    Unm1(kpos(1)), Unm1(kpos(2)), Kp);
                
    b_k = b_k + b_kloc;    

else % Covered: Kp lies completely in Gnm1_minus (xGl <= x_km1 & xGr <= b)
    
    b_k = [0; 0];
   
    % KGsCC, CC = Covered and Cut by Kp
    KGsCC = KGs(:, x_km1 < KGs(2,:) & KGs(1,:) < x_k);   
    for Km = KGsCC
        
        x_lm1 = Km(1);
        x_l = Km(2);
        x_lpos = find((G == x_l));
        lpos = I0 + (x_lpos - 1: x_lpos);
        
        x_le = max(x_lm1, x_km1); % Modify left endpoint if needed
        x_re = min(x_l, x_k); % Modify right endpoint if needed
        
        b_kloc = TimeJumpVector_Auxiliary_LocalofLocal_DiffMesh(x_le, ...
                    x_re, Unm1(lpos(1)), Unm1(lpos(2)), Kp, Km);
        
        b_k = b_k + b_kloc;
        
    end
    
end



