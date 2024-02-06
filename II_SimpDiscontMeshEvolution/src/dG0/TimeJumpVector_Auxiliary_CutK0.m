function b_k = TimeJumpVector_Auxiliary_CutK0(a, b, Kp, kpos, Unm1, KGs, G)

global I0

% Checks which simplices in mesh_composition at tnm1_minus that 
% potentially cut simplex Kp intersects and computes the local contribution 
% to rhs-vector by using integration interval (a, b), 
% where x_km1 <= a < b <= x_k.

% Note that this is m-file is a generalized and somewhat more expensive 
% version of the m-file "TimeJumpVector_Auxiliary_EntireK0".

% Simplex Kp belongs to mesh composition at time tnm1_plus
x_km1 = Kp(1);
x_k = Kp(2); 

% The left and right endpoints xGl and xGr belong to G(tnm1_minus)
xGl = G(1);
xGr = G(end);

% Check intersection and compute the contribution
if b <= xGl || xGr <= a % Uncovered: Potentially cut Kp is unchanged

    b_k = TimeJumpVector_Auxiliary_LocalofLocal_SameMesh(a, b, ...
                    Unm1(kpos(1)), Unm1(kpos(2)), Kp);
                   
elseif a < xGl && xGl < b % Cut: Potentially cut Kp contains xGl
 
    b_k = [0; 0];
    
    % Omega_1-part (interval: (a, xGl))
    b_kloc = TimeJumpVector_Auxiliary_LocalofLocal_SameMesh(a, xGl, ...
                    Unm1(kpos(1)), Unm1(kpos(2)), Kp);
                
    b_k = b_k + b_kloc;            
    
    % Omega_2/G-part (interval: (xGl, b))
    % KGsCC, CC = Covered and Cut by (a, b)
    KGsCC = KGs(:, KGs(1,:) < b);
    for Km = KGsCC
        
        x_lm1 = Km(1);
        x_l = Km(2);
        x_lpos = find((G == x_l));
        lpos = I0 + (x_lpos - 1: x_lpos);
        
        x_le = x_lm1;
        x_re = min(x_l, b); % Modify right endpoint if needed
        
        b_kloc = TimeJumpVector_Auxiliary_LocalofLocal_DiffMesh(x_le, ...
                    x_re, Unm1(lpos(1)), Unm1(lpos(2)), Kp, Km);
        
        b_k = b_k + b_kloc;
               
    end
        
elseif a < xGr && xGr < b % Cut: Potentially cut Kp contains xGr
 
    b_k = [0; 0];
    
    % Omega_2/G-part (interval: (a, xGr))
    % KGsCC, CC = Covered and Cut by (a, b)
    KGsCC = KGs(:, KGs(2,:) > a);
    for Km = KGsCC
      
        x_lm1 = Km(1);
        x_l = Km(2);
        x_lpos = find((G == x_l));
        lpos = I0 + (x_lpos - 1: x_lpos);
        
        x_le = max(x_lm1, a); % Modify left endpoint if needed
        x_re = x_l;
        
        b_kloc = TimeJumpVector_Auxiliary_LocalofLocal_DiffMesh(x_le, ...
                    x_re, Unm1(lpos(1)), Unm1(lpos(2)), Kp, Km);
        
        b_k = b_k + b_kloc;
        
    end
    
    % Omega_1-part (interval: (xGr, b))
    b_kloc = TimeJumpVector_Auxiliary_LocalofLocal_SameMesh(xGr, b, ...
                    Unm1(kpos(1)), Unm1(kpos(2)), Kp);
                
    b_k = b_k + b_kloc;    
    
else % Covered: Potentially cut Kp lies completely in Gnm1_minus 
     % (xGl<=a & b<=xGr)
 
    b_k = [0; 0];
   
    % KGsCC, CC = Covered and Cut by (a, b)
    KGsCC = KGs(:, a < KGs(2,:) & KGs(1,:) < b);   
    for Km = KGsCC
      
        x_lm1 = Km(1);
        x_l = Km(2);
        x_lpos = find((G == x_l));
        lpos = I0 + (x_lpos - 1: x_lpos);
        
        x_le = max(x_lm1, a); % Modify left endpoint if needed
        x_re = min(x_l, b); % Modify right endpoint if needed
        
        b_kloc = TimeJumpVector_Auxiliary_LocalofLocal_DiffMesh(x_le, ...
                    x_re, Unm1(lpos(1)), Unm1(lpos(2)), Kp, Km);
        
        b_k = b_k + b_kloc;
        
    end
    
end



