function b_fint = fIntegralTermdG1_BE(K0, KG, x0, G)

% Quadrature is used to approximate the integrals locally over 
% space-time prisms.
% 1st, temporal quadrature: 3-point Lobatto rule => quad error ~ k^4.
% 2nd, spatial quadrature: 3-point Gauss-Legendre rule => quad error ~ h^6.
% Old versions: 
% 2nd, spatial quadrature: Trapezoidal rule => quad error ~ h^2

global kn x0_init x0_fin I0 tnm1 tn M leA

b_fint = zeros(leA, 1);
bk = zeros(2, 1);

xGl = G(1);
xGr = G(end);

% Temporal quadrature points and "combined weights" = 3*qw_i*lamb_l(t_i)

% 3-point Lobatto quadrature in time
tm = 0.5*(tnm1 + tn);
t_vec = [tnm1, tm, tn];

Crs = [1 2 0; 0 2 1];

% % 3-point Gauss-Legendre quadrature in time
% tm = 0.5*(tnm1 + tn);
% t_vec = [tm - sqrt(3/5)*kn*0.5, tm, tm + sqrt(3/5)*kn*0.5];
% 
% Crs = (1/6)*[5+sqrt(15) 8 5-sqrt(15); 5-sqrt(15) 8 5+sqrt(15)];

% Domains of U1 ---------------------------------------------------

% To the left of G
K1ra = K0(:, K0(1,:) < xGl); 
leK1ra = length(K1ra(1,:));
for k = 1:leK1ra
    
    x_km1 = K1ra(1,k);
    x_k = K1ra(2,k);
    hk = K1ra(3, k);
    
    x_kpos = find((x0 == x_k));
    kpos = (x_kpos - 1: x_kpos);
    
    alpha = x_km1;
    beta = min(x_k, xGl);
    
    l = 1;
    bk1 = zeros(2,1);
    bk2 = zeros(2,1);
    for tl = t_vec
        
        % 3-point Gauss-Legendre rule in space
        bk(1) = Quad3pGL(@(x)(f_func(x, tl)*(x_k - x)), alpha, beta);
        bk(2) = Quad3pGL(@(x)(f_func(x, tl)*(x - x_km1)), alpha, beta);
        
%         % Trapezoidal rule in space
%         F = f_func(alpha,tl) + f_func(beta,tl); 
%         Fx = f_func(alpha,tl)*alpha + f_func(beta,tl)*beta;    
%         bk = [x_k*F - Fx; Fx - x_km1*F];    
    
        bk1 = bk1 + Crs(1,l)*bk;
        bk2 = bk2 + Crs(2,l)*bk;
        
        l = l + 1;
        
    end
       
    fac = (kn)/(6*hk); % 3-point Gauss-Legendre rule in space
%     fac = (beta - alpha)*kn/(12*hk); % Trapezoidal rule in space
    
    bk1 = fac*bk1;
    bk2 = fac*bk2;
    
    b_fint(kpos) = b_fint(kpos) + bk1;          % Test with lamb_nm1*phis
    b_fint(M + kpos) = b_fint(M + kpos) + bk2;  % Test with lamb_n*phis
    
end

% To the right of G
K1rb = K0(:, xGr < K0(2,:));
leK1rb = length(K1rb(1,:));
for k = 1:leK1rb
    
    x_km1 = K1rb(1,k);
    x_k = K1rb(2,k);
    hk = K1rb(3, k);
    
    x_kpos = find((x0 == x_k));
    kpos = (x_kpos - 1: x_kpos);
    
    alpha = max(x_km1, xGr);
    beta = x_k;
    
    l = 1;
    bk1 = zeros(2,1);
    bk2 = zeros(2,1);
    for tl = t_vec
        
        % 3-point Gauss-Legendre rule in space
        bk(1) = Quad3pGL(@(x)(f_func(x, tl)*(x_k - x)), alpha, beta);
        bk(2) = Quad3pGL(@(x)(f_func(x, tl)*(x - x_km1)), alpha, beta);
        
%         % Trapezoidal rule in space
%         F = f_func(alpha,tl) + f_func(beta,tl); 
%         Fx = f_func(alpha,tl)*alpha + f_func(beta,tl)*beta;    
%         bk = [x_k*F - Fx; Fx - x_km1*F];
    
        bk1 = bk1 + Crs(1,l)*bk;
        bk2 = bk2 + Crs(2,l)*bk;
        
        l = l + 1;
        
    end
    
    fac = (kn)/(6*hk); % 3-point Gauss-Legendre rule in space
%     fac = (beta - alpha)*kn/(12*hk); % Trapezoidal rule in space
    
    bk1 = fac*bk1;
    bk2 = fac*bk2;
    
    b_fint(kpos) = b_fint(kpos) + bk1;
    b_fint(M + kpos) = b_fint(M + kpos) + bk2;
    
end

% Domains of U2 -----------------------------------------------------------

K2r = KG(:, x0_init <= KG(1,:) & KG(2,:) <= x0_fin);
leK2r = length(K2r(1,:));
for k = 1:leK2r
    
    x_km1 = K2r(1,k);
    x_k = K2r(2,k);
    hk = K2r(3,k);
    
    x_kpos = find((G == x_k));
    kpos = I0 + (x_kpos - 1: x_kpos);
    
    alpha = x_km1;
    beta = x_k;
    
    l = 1;
    bk1 = zeros(2,1);
    bk2 = zeros(2,1);
    for tl = t_vec
        
        % 3-point Gauss-Legendre rule in space
        bk(1) = Quad3pGL(@(x)(f_func(x, tl)*(x_k - x)), alpha, beta);
        bk(2) = Quad3pGL(@(x)(f_func(x, tl)*(x - x_km1)), alpha, beta);
        
%         % Trapezoidal rule in space
%         F = f_func(alpha,tl) + f_func(beta,tl); 
%         Fx = f_func(alpha,tl)*alpha + f_func(beta,tl)*beta;    
%         bk = [x_k*F - Fx; Fx - x_km1*F];        
    
        bk1 = bk1 + Crs(1,l)*bk;
        bk2 = bk2 + Crs(2,l)*bk;
        
        l = l + 1;
        
    end
    
    fac = (kn)/(6*hk); % 3-point Gauss-Legendre rule in space
%     fac = kn/12; % Trapezoidal rule in space
    
   bk1 = fac*bk1;
   bk2 = fac*bk2;
    
   b_fint(kpos) = b_fint(kpos) + bk1;
   b_fint(M + kpos) = b_fint(M + kpos) + bk2;
    
end  