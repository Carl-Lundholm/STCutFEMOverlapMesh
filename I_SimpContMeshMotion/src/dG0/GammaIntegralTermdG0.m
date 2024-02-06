function A_git = GammaIntegralTermdG0(K0, KG, x0, mutr, t1ga, t1gb, ...
                                amin, amax, bmin, bmax, w1, w2, pnltp)

global  I0 IG leA

C = sqrt(1 + mutr^2);
D = 1/C;

A_git = zeros(leA);

% Compute the cut the K0:s
if mutr == 0 
    % Special conditions for left and right of G when mu = 0
    Kgs_left = K0(:, K0(1,:) < amin & amax <= K0(2,:));
    Kgs_right = K0(:, K0(1,:) <= bmin & bmax < K0(2,:));
else 
    % General conditions for left and right of G
    Kgs_left = K0(:, K0(1,:) < amax & amin < K0(2,:));
    Kgs_right = K0(:, K0(1,:) < bmax & bmin < K0(2,:));
end

% To the left of G --------------------------------------------------------
n1 = 1;
nt = n1*mutr*D;     % The extra minus sign from the FE-form is included 
nx = -n1*D;         % - || -

GbL = I0 + 1;
GbR = I0 + 2;
Gb = GbL;
hGb = KG(3,1);

tg = t1ga;
         
A_git = GaInTedG0_AuxiliaryFunc(A_git, Kgs_left, x0, C, tg, nt, nx, ...
    GbL, GbR, Gb, hGb, amin, amax, w1, w2, pnltp);
  
% To the right of G -------------------------------------------------------
n1 = -1;
nt = n1*mutr*D;     % The extra minus sign from the FE-form is included 
nx = -n1*D;         % - || -

GbL = leA - 1;
GbR = leA;
Gb = GbR;
hGb = KG(3, IG - 1);

tg = t1gb;

A_git = GaInTedG0_AuxiliaryFunc(A_git, Kgs_right, x0, C, tg, nt, nx, ...
    GbL, GbR, Gb, hGb, bmin, bmax, w1, w2, pnltp);
  
       
 
    
    
    
    
    
    
    
    
    
    
    