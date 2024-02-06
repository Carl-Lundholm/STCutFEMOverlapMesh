function b_fint = fIntegralTermdG0_BE(K0, KG, x0, G)

% Quadrature is used to approximate the integrals locally over 
% space-time prisms.
% 1st, temporal quadrature: midpoint rule => quad error ~ k^2.
% 2nd, spatial quadrature: 3-point Gauss-Legendre rule => quad error ~ h^6.
% Old versions: 
% 2nd, spatial quadrature: Trapezoidal rule => quad error ~ h^2

global kn x0_init x0_fin I0 tnm1 tn IG

b_fint = zeros(I0 + IG, 1);
bk = zeros(2, 1);

xGl = G(1);
xGr = G(end);

t_av = 0.5*(tn + tnm1);         % Time that is used in the quadrature

% Domains of U1 ---------------------------------------------------

% To the left of G
K1ra = K0(:, K0(1,:) < xGl); 
leK1ra = length(K1ra(1,:));
for k = 1:leK1ra
    
    x_km1 = K1ra(1,k);
    x_k = K1ra(2,k);
    hk = K1ra(3,k);
    
    x_kpos = find((x0 == x_k));
    kpos = (x_kpos - 1: x_kpos);
    
    alpha = x_km1;
    beta = min(x_k, xGl);
    
    % 3-point Gauss-Legendre rule in space
    bk(1) = Quad3pGL(@(x)(f_func(x, t_av)*(x_k - x)), alpha, beta);
    bk(2) = Quad3pGL(@(x)(f_func(x, t_av)*(x - x_km1)), alpha, beta);
    
    b_fint(kpos) = b_fint(kpos) + (kn/hk)*bk;
    
%     % Trapezoidal rule in space
%     delx = beta - alpha;
%      
%     alpha_kn = 0.5*(f_func(alpha,t_av) + f_func(beta,t_av));
%     beta_kn = 0.5*(f_func(alpha,t_av)*alpha + f_func(beta,t_av)*beta);
%     
%     bk = 1/hk*delx*kn*[x_k*alpha_kn - beta_kn; beta_kn - x_km1*alpha_kn];
% 
%     b_fint(kpos) = b_fint(kpos) + bk;
    
end

% To the right of G
K1rb = K0(:, xGr < K0(2,:));
leK1rb = length(K1rb(1,:));
for k = 1:leK1rb
    
    x_km1 = K1rb(1,k);
    x_k = K1rb(2,k);
    hk = K1rb(3,k);
    
    x_kpos = find((x0 == x_k));
    kpos = (x_kpos - 1: x_kpos);
    
    alpha = max(xGr, x_km1);
    beta = x_k;
    
    % 3-point Gauss-Legendre rule in space
    bk(1) = Quad3pGL(@(x)(f_func(x, t_av)*(x_k - x)), alpha, beta);
    bk(2) = Quad3pGL(@(x)(f_func(x, t_av)*(x - x_km1)), alpha, beta);
    
    b_fint(kpos) = b_fint(kpos) + (kn/hk)*bk;
    
%     % Trapezoidal rule in space
%     delx = beta - alpha;
%      
%     alpha_kn = 0.5*(f_func(alpha,t_av) + f_func(beta,t_av));
%     beta_kn = 0.5*(f_func(alpha,t_av)*alpha + f_func(beta,t_av)*beta);
%     
%     bk = 1/hk*delx*kn*[x_k*alpha_kn - beta_kn; beta_kn - x_km1*alpha_kn];
% 
%     b_fint(kpos) = b_fint(kpos) + bk;
    
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
    
    % 3-point Gauss-Legendre rule in space
    bk(1) = Quad3pGL(@(x)(f_func(x, t_av)*(x_k - x)), x_km1, x_k);
    bk(2) = Quad3pGL(@(x)(f_func(x, t_av)*(x - x_km1)), x_km1, x_k);
    
    b_fint(kpos) = b_fint(kpos) + (kn/hk)*bk;
    
%     % Trapezoidal rule in space
%     alpha_kn = 0.5*(f_func(x_km1,t_av) + f_func(x_k,t_av));
%     beta_kn = 0.5*(f_func(x_km1,t_av)*x_km1 + f_func(x_k,t_av)*x_k);
%     
%     bk = kn*[x_k*alpha_kn - beta_kn; beta_kn - x_km1*alpha_kn];
%     
%     b_fint(kpos) = b_fint(kpos) + bk;
    
end
    
    
    
    
    
    
    
    
    
    