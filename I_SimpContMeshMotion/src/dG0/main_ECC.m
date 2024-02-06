% clc
% clear all
% close fig

function [error, kn, h0] = main_ECC(kh_init, hkfix)

global I0 IG I1 I2 x0_init x0_fin t0 a b h0 hG kn anm1 bnm1 ...
       kn1ga kn1gb tnm1 tn mutr leA G Gnm1 M N x0 K0 KG ...
       amin amax bmin bmax t1ga t1gb C w1 w2

k_ECC = kh_init;
h_ECC = k_ECC;

% parameters_ECC
parameters_ECC_k_fixedh % Use k_ECC to compute kn and use fixed h
% parameters_ECC_h_fixedk % Overwrite kn with fixed number and use h = h_ECC = k_ECC

x0 = (x0_init:h0:x0_fin);
K0 = [x0(1:I0-1); x0(2:I0); x0(2:I0) - x0(1:I0-1)];

G = (a:hG:b);
KG = [G(1:IG-1); G(2:IG); G(2:IG) - G(1:IG-1)];

tn = t0;
delxh0 = kn*mutr/h0;

disp(['mu = ', num2str(mutr)]);
disp(['Number of time steps = ', num2str(N - 1)]);
disp(['h_ECC = ', num2str(h_ECC)]);
disp(['hG = ', num2str(hG)]);
disp(['h0 = ', num2str(h0)]);
disp(['kn = ', num2str(kn)]);
% disp(['delx = ', num2str(kn*mutr)]);
% disp(['pnltp = ', num2str(pnltp)]);
% disp(' ');
% 
% disp('Convergence preserving criteria');
% 
% disp('abs(mu) < h/kn ?');
% h = max(hG,h0);
% if abs(mutr) < h/kn
%     disp('yes');
% else 
%     disp('no');
% end
% 
% disp('h^2 <= kn ?');
% if h^2 <= kn
%     disp('yes');
% else 
%     disp('no');
% end
% disp(' ');

anm1 = a;
bnm1 = b;
[x1, x2, infoOmega_2] = MeshGenerator(x0, G);

I1 = length(x1);
I2 = length(x2);

U0a = u_func(x0,tn);
U2a = u_func(G,tn);
U = [U0a, U2a]';

U_q2p = U;

M = I0 + IG;
leA = M;

error_Xnorm = 0;

for n=2:N
    
    tnm1 = tn;
    tn = tn + kn;
    
    Gnm1 = G;
    KGnm1 = KG;
    anm1 = a;
    bnm1 = b;
    
    mutr = mu_func(tn, mutr, a, b);
    C = sqrt(1 + mutr^2);
    G = Gnm1 + mutr*kn;
    KG = [G(1:IG-1); G(2:IG); G(2:IG) - G(1:IG-1)];
      
    x1nm1 = x1;
    x2nm1 = x2;
    
    Unm1 = U;
    Unm1_q2p = U_q2p;
    
    a = G(1);
    b = G(IG);
    
    [x1, x2, infoOmega_2] = MeshGenerator(x0, G);                             
    
    amin = min(anm1,a);
    amax = max(anm1,a);
    bmin = min(bnm1,b);
    bmax = max(bnm1,b);
    
    x1r = x1(x1 <= amin | bmax <= x1);
    
    x1ga = x1(amin < x1 & x1 < amax);
    x1gb = x1(bmin < x1 & x1 < bmax);
    
    t1ga = tnm1 + (x1ga - anm1)/mutr;
    t1gb = tnm1 + (x1gb - bnm1)/mutr; 
    
    if mutr >= 0
        
        t1ga = [tnm1, t1ga, tn];
        t1gb = [tnm1, t1gb, tn];
        
        tna = tn*ones(1,length(t1ga));
        tnm1a = t1ga;
        tnb = t1gb;
        tnm1b = tnm1*ones(1,length(t1gb));
    else
        
        t1ga = [tn, t1ga, tnm1];
        t1gb = [tn, t1gb, tnm1];
        
        tna = t1ga;
        tnm1a = tnm1*ones(1,length(t1ga));
        tnb = tn*ones(1,length(t1gb));
        tnm1b = t1gb;
    end
    
    kn1ga = [tna - tnm1a; 0.5*(tna + tnm1a)];
    kn1gb = [tnb - tnm1b; 0.5*(tnb + tnm1b)];  
                
    I1 = length(x1);
    I2 = length(x2);
    
    % Analytic solution ---------------------------------------------------   

    U0a = u_func(x0,tn);
    U2a = u_func(x2,tn);
    Ua = [U0a, U2a]';
             
    % Computation of A_mat ------------------------------------------------
    
    % Initialize matrix A_mat
    A_mat = zeros(leA);
    
    % Contribution to A_mat from the ALE term
    A_ale = ALEdG0(mutr, kn);
    
    A_mat = A_mat + A_ale;
    
    % Contribution to A_mat from the time-jump terms
    A_tjnm1 = TimeJumpMatrix(K0, KG, x0, G, anm1, bnm1);  
    
    A_mat = A_mat + A_tjnm1;
    
    % Contribution to A_mat from the stiffness matrix     
%     [A_sm, b_fint] = StiffnessMatrixdG0fIntegral(K0, KG, x0, G, ...
%                               amin, amax, bmin, bmax);
    A_sm = StiffnessMatrixdG0(K0, KG, x0, G, ...
                                amin, amax, bmin, bmax);
    %A_sm_q2p = StiffnessMatrixdG0_q2p(K0, KG, x0, G, ...
    %                            amin, amax, bmin, bmax);
    
    A_mat = A_mat + A_sm;

    % Contribution to A_mat from the gamma integral terms
    A_git = GammaIntegralTermdG0(K0, KG, x0, mutr, t1ga, t1gb, ...
                               amin, amax, bmin, bmax, w1, w2, pnltp);
                           
    A_mat = A_mat + A_git;
    
    % Contribution to A_mat from the derivative penalty term
    A_pnt = PenaltyNablaTermdG0(K0, KG, x0, mutr, t1ga, t1gb, ...
                                amin, amax, bmin, bmax, pnltp_nab);
                            
    A_mat = A_mat + A_pnt;            
    
    % Contribution to A_mat from the boundary conditions -----
%     % Add Robin boundary conditions 
%     A_bcx0init = -eta*kn;
%     A_bcx0fin = A_bcx0init;
%     
%     A_mat(1, 1) = A_mat(1, 1) + A_bcx0init;
%     A_mat(I0, I0) = A_mat(I0, I0) + A_bcx0fin;
    
    % Apply Dirichlet boundary conditions   
    A_mat(1, :) = zeros(leA, 1);
    A_mat(1, 1) = 1;

    A_mat(I0, :) = zeros(leA, 1);
    A_mat(I0, I0) = 1;
    
    % Computation of b_vec ------------------------------------------------
    
    % Initialize vector b
    b_vec = zeros(leA, 1);
    
    % Contribution to b_vec from the time-jump terms
    b_tjnm1 = A_tjnm1*Unm1;
    %b_tjnm1_q2p = A_tjnm1*Unm1_q2p;
    
    b_vec = b_vec + b_tjnm1;  
    
    % Contribution to b_vec from the interior f-integral terms
    b_fint = fIntegralTermdG0(K0, KG, x0, G, ...
                               amin, amax, bmin, bmax, mutr, t1ga, t1gb);
%     b_fint = fIntegralTermdG0_cruderQuad(K0, KG, x0, G, ...
%                                 amin, amax, bmin, bmax);
    b_vec = b_vec + b_fint;
    
    % Contribution to b_vec from the boundary conditions -----
%     % Add Robin boundary conditions
%     Integrand_bc = @(x,t)(gN_func(x,t) - eta*gD_func(x,t));
%     % Computation with MATLAB's INTEGRAL    
% %     b_bcx0init = integral(@(t)Integrand_bc(x0_init, t), tnm1, tn);
% %     b_bcx0fin = integral(@(t)Integrand_bc(x0_fin, t), tnm1, tn);
% 
%     % Computation with quadrature
%     b_bcx0init = Integrand_bc(x0_init, (tn + tnm1)/2)*kn;
%     b_bcx0fin = Integrand_bc(x0_fin, (tn + tnm1)/2)*kn;
% 
%     b_vec(1) = b_vec(1) + b_bcx0init;
%     b_vec(I0) = b_vec(I0) + b_bcx0fin;
    
    % Apply Dirichlet boundary conditions    
    Ubc = 0;      % u on the boundary of Omega_0
    
    b_vec(1) = Ubc;  
    b_vec(I0) = Ubc;   
    
    % Solving the linear equation system, A_mat*U = b_vec -----------------
    
    acno0 = ones(leA,1);
    acno = acno0;
    for k = 1:I0
        
        if ~any(A_mat(k,:)) 
           A_mat(k, k) = 1;
           b_vec(k) = 0;
           Ua(k) = 0;
           acno(k) = 0;
        end
        
    end

    U = A_mat\b_vec; 
    
    % Solving the linear equation system, A_mat_q2p*U_q2p = b_vec_q2p ----- 
   
    %A_mat_q2p = A_mat - A_sm + A_sm_q2p;
    %b_vec_q2p = b_vec - b_tjnm1 + b_tjnm1_q2p;
    
    %U_q2p = A_mat_q2p\b_vec_q2p;
    
    % Compute the contribution to the squared X-norm, i.e., space-time
    % energy norm, of the error from the current space-time slab
    error_Xnorm_slab = ErrorXNormSlabwiseSquared(U, U, Unm1, @u_func, n);
                
    error_Xnorm = error_Xnorm + error_Xnorm_slab;    

end
    
% Computation of errors ---------------------------------------------------

% X-norm of error: Take the square root
error_Xnorm = sqrt(error_Xnorm);

% Return the X-norm of the error
error = error_Xnorm;

% Return the L2-norm of the error at the final time
% errexa = ErrorL2NormSpaceQuad3PointGaussLegendre(U, ...
%                                     @(x)(u_func(x,tn)), x0, G);
% 
% %errexa = ErrorL2NormSpaceQuadTrapezoidal(U, Ua, x0, G, a, b);
% 
% ex_err = Ua(acno0 == acno) - U(acno0 == acno);
% errexa_approx = sqrt(mean(abs(ex_err).^2)*(x0_fin - x0_init)); 
% 
% % Quadrature q2p err, Ua - U_q2p
% errL2Norm_q2p = ErrorL2NormSpaceQuadTrapezoidal(U_q2p, Ua, x0, G, a, b);
% errq2p = errL2Norm_q2p;
% 
% q2p_err = Ua(acno0 == acno) - U_q2p(acno0 == acno);
% errq2p_approx = sqrt(mean(abs(q2p_err).^2)*(x0_fin - x0_init));    








    
    
    