% clc
% clear all
% close fig

function [error, kn, h0] = MMNP_main_ECC(kh_init, mu_iter)

global I0 IG I1 I2 x0_init x0_fin t0 a b h0 hG kn anm1 bnm1 ...
       kn1ga kn1gb tnm1 tn mutr leA G Gnm1 M Ind N x0 K0 KG ...
       amin amax bmin bmax t1ga t1gb C w1 w2

k_ECC = kh_init;
h_ECC = k_ECC;

MMNP_parameters_ECC_k_fixedh % Use k_ECC to compute kn and use fixed h
% MMNP_parameters_ECC_h_fixedk % Overwrite kn with fixed number and use h = h_ECC = k_ECC

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
% disp('kn*abs(mu) < h ?');
% h = max(hG,h0);
% if kn*abs(mutr) < h
%     disp('yes');
% else 
%     disp('no');
% end
% 
% disp('h^1 <= kn^3/2 ?');
% if h^1 <= kn^(3/2)
%     disp('yes');
% else 
%     disp('no');
% end
% disp(' ');
% 
% disp('h^2 <= kn^3 ?');
% if h^2 <= kn^3
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

U0an = u_func(x0,tn);
U2an = u_func(G,tn);
Uan = [U0an, U2an]';
Un = Uan;

U_q2p = Un;

error_Xnorm = 0;

M = I0 + IG;
leA = 2*M; 
Ind = [0 M];

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
    
    Unm1 = Un;
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
    
    kn1ga = [tna - tnm1a; (tna + tnm1a)/2];
    kn1gb = [tnb - tnm1b; (tnb + tnm1b)/2];  
                
    I1 = length(x1);
    I2 = length(x2);
    
    % Analytic solution ---------------------------------------------------   
    Uanm1 = Uan;
    
    U0an = u_func(x0,tn);
    U2an = u_func(G,tn);
    Uan = [U0an, U2an]';
             
    % Computation of A_mat ------------------------------------------------
    
    % Initialize matrix A_mat
    A_mat = zeros(leA);
    
    % Contribution to A_mat from the D_t-term -----
    A_dt = DtTermdG1(K0, KG, x0, G, ...
                                amin, amax, bmin, bmax, t1ga, t1gb);
    
    A_mat = A_mat + A_dt;
    
    % Contribution to A_mat from the ALE term -----
    A_ale = ALEdG1(KG, G, mutr, kn);
    
    A_mat = A_mat + A_ale;
    
    % Contribution to A_mat from the time jump term -----
    A_tjnm1 = TimeJumpMatrix(K0, KG, x0, G, anm1, bnm1);  
    
    A_mat(1:M,1:M) = A_mat(1:M,1:M) + A_tjnm1;
    
    % Contribution to A_mat from the stiffness matrix -----     
%     [A_sm, b_fint] = StiffnessMatrixdG1fIntegral(K0, KG, x0, G, ...
%                               amin, amax, bmin, bmax);
    A_sm = StiffnessMatrixdG1(K0, KG, x0, G, ...
                            amin, amax, bmin, bmax, t1ga, t1gb);
%     A_sm_q2p = StiffnessMatrixdG1_q2p(K0, KG, x0, G, ...
%         amin, amax, bmin, bmax);

    A_mat = A_mat + A_sm;

    % Contribution to A_mat from the gamma integral terms -----
    A_git = GammaIntegralTermdG1(K0, KG, x0, mutr, C, t1ga, t1gb, ...
                               amin, amax, bmin, bmax, w1, w2, pnltp);
    
    A_mat = A_mat + A_git;
    
    % Contribution to A_mat from the derivative penalty term -----
    A_pnt = PenaltyNablaTermdG1(K0, KG, x0, mutr, C, t1ga, t1gb, ...
                                amin, amax, bmin, bmax, pnltp_nab);
                            
    A_mat = A_mat + A_pnt;            
    
    % Contribution to A_mat from the boundary conditions  -----   
%     % Add Robin boundary conditions 
%     A_bc = -eta*kn/6;
%     
%     Crs = [2 1; 1 2];
%  
%     A_mat(Ind + 1, Ind + 1) = A_mat(Ind + 1, Ind + 1) + A_bc*Crs;
%     A_mat(Ind + I0, Ind + I0) = A_mat(Ind + I0, Ind + I0) + A_bc*Crs;

    % Apply Dirichlet boundary conditions   
    for r = Ind
        A_mat(r + 1, :) = zeros(leA, 1);
        A_mat(r + 1, r + 1) = 1;
        
        A_mat(r + I0, :) = zeros(leA, 1);
        A_mat(r + I0, r + I0) = 1;
    end              
    
    % Computation of b_vec ------------------------------------------------
    
    % Initialize vector b
    b_vec = zeros(leA, 1);
    
    % Contribution to b_vec from the time-jump term -----
    b_tjnm1 = A_tjnm1*Unm1;
    % b_tjnm1_q2p = A_tjnm1*Unm1_q2p;
    
    b_vec(1:M) = b_vec(1:M) + b_tjnm1;  
    
    % Contribution to b_vec from the interior f-integral terms -----
    b_fint = fIntegralTermdG1(K0, KG, x0, G, ...
                            amin, amax, bmin, bmax, mutr, t1ga, t1gb);
    
    b_vec = b_vec + b_fint;
    
    % Contribution to b_vec from the boundary conditions -----   
%     % Add Robin boundary conditions 
%     Integrand_bc = @(x,t)(gN_func(x,t) - eta*gD_func(x,t));
%     % Computation with MATLAB's INTEGRAL    
% %     b_bcx0init = integral(@(t)Integrand_bc(x0_init, t), tnm1, tn);
% %     b_bcx0fin = integral(@(t)Integrand_bc(x0_fin, t), tnm1, tn);
% 
%     % Computation with 3 point Lobatto quadrature
%     tm = 0.5*(tnm1 + tn);
%     t_vec = [tnm1 tm tn];
%     
%     Crs = [1 2 0; 0 2 1];
%     l = 1;
%     
%     b_bcx0init = zeros(2,1);
%     b_bcx0fin = zeros(2,1);
%     for tl = t_vec;
%         
%         intbc_x0init_tl = Integrand_bc(x0_init, tl);
%         intbc_x0fin_tl = Integrand_bc(x0_fin, tl);
%         
%         b_bcx0init(1) = b_bcx0init(1) + Crs(1,l)*intbc_x0init_tl;
%         b_bcx0init(2) = b_bcx0init(2) + Crs(2,l)*intbc_x0init_tl;
%         
%         b_bcx0fin(1) = b_bcx0fin(1) + Crs(1,l)*intbc_x0fin_tl;
%         b_bcx0fin(2) = b_bcx0fin(2) + Crs(2,l)*intbc_x0fin_tl;
%     
%         l = l + 1;
%     end
%    
%     b_vec(Ind + 1) = b_vec(Ind + 1) + kn/6*b_bcx0init;
%     b_vec(Ind + I0) = b_vec(Ind + I0) + kn/6*b_bcx0fin; 
    
    % Apply Dirichlet boundary conditions    
    Ubc = 0;      % u on the boundary of Omega_0
    
    for r = Ind
        b_vec(r + 1) = Ubc;  
        b_vec(r + I0) = Ubc;
    end
    
    % Solving the linear equation system, A_mat*U = b_vec -----------------
       
    % Fixing the natural gaps along the diagonal of A_mat where G lies
    acno0 = ones(leA,1);
    acno = acno0; 
    for r = Ind 
        for k = 1:I0
            if ~any(A_mat(r + k, r + 1: r + M)) 
                A_mat(r + k, r + k) = 1;
                b_vec(r + k) = 0;
                if r == M
                    Uan(k) = 0;
                end
                acno(r + k) = 0;
            end
        end
    end
    
    U = A_mat\b_vec; 
    
    Unm1p = U(1:M);
    Un = U(M+1:leA);
    
    % Compute the contribution to the squared X-norm, i.e., space-time
    % energy norm, of the error from the current space-time slab
    error_Xnorm_slab = ErrorXNormSlabwiseSquared(Un, Unm1p, Unm1, ...
                            @u_func, n);
                
    error_Xnorm = error_Xnorm + error_Xnorm_slab;             

end
    
% Computation of errors ---------------------------------------------------
% X-norm of error: Take the square root
error_Xnorm = sqrt(error_Xnorm);

% Return the X-norm of the error
error = error_Xnorm;
  
% Return the L2-norm of the error the the final time
%error = ErrorL2NormSpaceQuadTrapezoidal(Un, Uan, x0, G, a, b);    
%error = ErrorL2NormSpaceQuad(Un, @(x)(u_func(x, tn)), x0, G);
    







    
    
    