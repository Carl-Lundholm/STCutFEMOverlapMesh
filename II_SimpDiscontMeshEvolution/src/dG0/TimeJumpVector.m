function b_tjnm1 = TimeJumpVector(K0, KG, x0, G, Unm1, KGnm1, Gnm1)

global x0_init x0_fin I0 IG

b_tjnm1 = zeros(I0 + IG, 1);

xGl = G(1);
xGr = G(end);

% Regular domains of U1 ---------------------------------------------------

% Entire/Uncut simplices K0 (regular domains to the left and right of G)
K0s_entire = K0(:, K0(2,:) <= xGl | xGr <= K0(1,:));
for Kp = K0s_entire
   
    x_k = Kp(2);
    x_kpos = find((x0 == x_k));
    kpos = (x_kpos - 1: x_kpos);
    
    b_k = TimeJumpVector_Auxiliary_EntireK0(Kp, kpos, ...
        Unm1, KGnm1, Gnm1);
    
    b_tjnm1(kpos, 1) = b_tjnm1(kpos, 1) + b_k;
    
end

% Boundary domains of U1 --------------------------------------------------

% To the left of G
Kp = K0(:, K0(1,:) < xGl & xGl < K0(2,:));
if ~isempty(Kp)
    
    x_km1 = Kp(1);
    x_k = Kp(2);
    x_kpos = find((x0 == x_k));
    kpos = (x_kpos - 1: x_kpos);

    b_k = TimeJumpVector_Auxiliary_CutK0(x_km1, xGl, Kp, ... 
        kpos, Unm1, KGnm1, Gnm1);
    
    b_tjnm1(kpos, 1) = b_tjnm1(kpos, 1) + b_k;
end


% To the right of G
Kp = K0(:, K0(1,:) < xGr & xGr < K0(2,:));
if ~isempty(Kp)
    
    x_km1 = Kp(1);
    x_k = Kp(2);
    x_kpos = find((x0 == x_k));
    kpos = (x_kpos - 1: x_kpos);

    b_k = TimeJumpVector_Auxiliary_CutK0(xGr, x_k, Kp, ...
        kpos, Unm1, KGnm1, Gnm1);
    
    b_tjnm1(kpos, 1) = b_tjnm1(kpos, 1) + b_k;
end


% Domains of U2 -----------------------------------------------------------
K2 = KG(:, x0_init <= KG(1,:) & KG(2,:) <= x0_fin);
for Kp = K2
   
    x_k = Kp(2);
    x_kpos = find((G == x_k));
    kpos = I0 + (x_kpos - 1: x_kpos);

    b_k = TimeJumpVector_Auxiliary_EntireKG(Kp, ...
        Unm1, K0, KGnm1, x0, Gnm1);
    
    b_tjnm1(kpos, 1) = b_tjnm1(kpos, 1) + b_k;  
    
end



