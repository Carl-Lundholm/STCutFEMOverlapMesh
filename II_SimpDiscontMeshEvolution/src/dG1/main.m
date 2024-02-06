
clc
clear all
close all

global I0 IG I1 I2 x0_init x0_fin t0 a b h0 hG kn ...
       kn1ga kn1gb tnm1 tn leA G M Ind

% parameters
% parameters_debug_dG1_ALogg
% parameters_debug_dG1
% parameters_0G_nodecollision
% parameters_manML
% parameters_manML_crude
% parameters_goodexofNabPen
% parameters_trapfornicecode
% parameters_nice_mu5em2
% parameters_PSACT
% parameters_boljfuncyta
% parameters_mesh_coarse
parameters_mesh_fine
% parameters_termtests

x0 = (x0_init:h0:x0_fin);
K0 = [x0(1:I0-1); x0(2:I0); x0(2:I0) - x0(1:I0-1)];

G = (a:hG:b);
KG = [G(1:IG-1); G(2:IG); G(2:IG) - G(1:IG-1)];
  
N = ceil((T - t0)/kn) + 1;
tn = t0;

disp('------------------------------------------------------------------');
disp(['Number of time steps = ', num2str(N - 1)]);
disp(['hG = ', num2str(hG)]);
disp(['h0 = ', num2str(h0)]);
disp(['delx = ', num2str(kn*mutr)]);
disp(['pnltp = ', num2str(pnltp)]);
disp(' ');

[x1, x2, infoOmega_2] = MeshGenerator(x0, G);

I1 = length(x1);
I2 = length(x2);

U0an = u_func(x0,tn);
U2an = u_func(G,tn);
Uan = [U0an, U2an]';
Un = Uan;

f0n = f_func(x0,tn);
f2n = f_func(x2,tn);
fn = [f0n, f2n]';

M = I0 + IG;
leA = 2*M; 
Ind = [0 M];

for n=2:N
  
    tnm1 = tn;
    tn = tn + kn;
    
    Gnm1 = G;
    KGnm1 = KG;
    Unm1 = Un;    
    
    mutr = mu_func(tn, mutr, a, b);
    C = sqrt(1 + mutr^2);
    G = G + mutr*kn;
    KG = [G(1:IG-1); G(2:IG); G(2:IG) - G(1:IG-1)];
      
    a = G(1);
    b = G(IG);
 
    [x1, x2, infoOmega_2] = MeshGenerator(x0, G);                             
    
    x1nm1 = x1;
    x2nm1 = x2;
    
    x1r = x1(x1 <= a | b <= x1);
    
    x1ga = [];
    x1gb = [];
    
    if mutr >= 0
        
        t1ga = [tnm1, tn];
        t1gb = [tnm1, tn];
        
        tna = tn*ones(1,length(t1ga));
        tnm1a = t1ga;
        tnb = t1gb;
        tnm1b = tnm1*ones(1,length(t1gb));
    else
        
        t1ga = [tn, tnm1];
        t1gb = [tn, tnm1];
        
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
    
    % The f function ------------------------------------------------------   
    fnm1 = fn;
    
    f0n = f_func(x0,tn);
    f2n = f_func(G,tn);
    fn = [f0n, f2n]';
    
    % Mesh Plot -----------------------------------------------------------
%     figure(1)
%     PlotMesh(G, x1, x2, x1r, x1ga, x1gb, x2nm1, ...
%              tn, tnm1, tna, tnm1a, tnb, tnm1b)
%     hold on
%     
%     axis([x0_init x0_fin t0 T])
%     xlabel('$x$', 'interpreter', 'latex', 'Fontsize', 24);
%     ylabel('$t$', 'interpreter', 'latex', 'Fontsize', 24);
    
%     figure(2)
%     PlotMesh(x1, x2, x1r, x1ga, x1gb, x2nm1, ...
%              tn, tnm1, tna, tnm1a, tnb, tnm1b)
%     hold on

%     figure(1);
%     plot(x1, zeros(I1,1), 'bo-');
%     hold on
%     if ~isempty(x2)
%         plot(x2, zeros(I2,1), 'rx-');
%     end
%     hold on
             
    % Computation of A_mat ------------------------------------------------
    
    % Initialize matrix A_mat
    A_mat = zeros(leA);
    
    % Contribution to A_mat from the delta_t term -----
    A_dt = DeltatTermdG1_BE(K0, KG, x0, G);
    
    A_mat = A_mat + A_dt;
   
    % BE-MESH-SETTING NOTE: There is no ALE-term for this setting. 
    
    % Contribution to A_mat from the time jump term -----
    A_tjnm1 = TimeJumpMatrix_BE(K0, KG, x0, G);  
    
    A_mat(1:M,1:M) = A_mat(1:M,1:M) + A_tjnm1;
    
    % Contribution to A_mat from the stiffness matrix -----     
    A_sm = StiffnessMatrixdG1_BE(K0, KG, x0, G);

    A_mat = A_mat + A_sm;

    % Contribution to A_mat from the gamma integral terms -----
    A_git = GammaIntegralTermdG1_BE(K0, KG, x0, t1ga, t1gb, ...
                               a, b, w1, w2, pnltp);
    
    A_mat = A_mat + A_git;
    
    % Contribution to A_mat from the derivative penalty term -----
    A_pnt = PenaltyNablaTermdG1_BE(K0, KG, x0, t1ga, t1gb, ...
                                a, b, pnltp_nab);
                            
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
    
    % BE-MESH-SETTING NOTE: Below is one of the biggest differences
    % compared to having a continuous mesh movement. With a continuous mesh
    % movement the time-jump matrix can be used here as well since the mesh
    % compositions at times tnm1_minus and tnm1_plus coincide. This is in
    % general not true for the BE-mesh-setting, hence the need for the
    % function TimeJumpVector that assembles (u_{h,n-1}^-, v_{n-1}^+).
    
    % Contribution to b_vec from the time-jump term -----
    b_tjnm1 = TimeJumpVector(K0, KG, x0, G, Unm1, KGnm1, Gnm1); 
    
    %b_tjnm1 = A_tjnm1*Unm1;
    % b_tjnm1_q2p = A_tjnm1*Unm1_q2p;
    
    b_vec(1:M) = b_vec(1:M) + b_tjnm1;  
    
    % Contribution to b_vec from the interior f-integral terms -----
    b_fint = fIntegralTermdG1_BE(K0, KG, x0, G);
    
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
   
%     cond(A_mat)
    
    U = A_mat\b_vec; 
    
    Unm1p = U(1:M);
    Un = U(M+1:leA);
    
    % Solution Plot -------------------------------------------------------
       
    figure(2)
    PlotSolutiondG1(Un, Unm1p, x0, G, G, a, a, b, b, ...
                                    tna, tnm1a, tnb, tnm1b)
    hold on                            
    axis([x0_init x0_fin t0 T])
    xlabel('$x$', 'interpreter', 'latex', 'Fontsize', 24);
    ylabel('$t$', 'interpreter', 'latex', 'Fontsize', 24);
    zlabel('$u_h$', 'interpreter', 'latex', 'Fontsize', 24)
    
end

% Add axis properties -----------------------------------------------------

axis([x0_init x0_fin t0 T])
xlabel('$x$', 'interpreter', 'latex', 'Fontsize', 24);
ylabel('$t$', 'interpreter', 'latex', 'Fontsize', 24);
zlabel('$u_h$', 'interpreter', 'latex', 'Fontsize', 24)
    
% Error computation -------------------------------------------------------

% FE-Sol err, Ua - U
% error = ErrorL2NormSpaceQuadTrapezoidal(Un, Uan, x0, G, a, b);
error = ErrorL2NormSpace(Un, @(x)(u_func(x,tn)), x0, G);

disp('L2-norm of the error of the dG(1)-FE-solution'); 
disp([num2str(error)]);
disp(' ');




    
    
    